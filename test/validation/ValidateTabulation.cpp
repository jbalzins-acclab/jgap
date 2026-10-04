#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include <highfive/H5File.hpp>
#include <highfive/H5Group.hpp>

#include "common/ValidationUtils.hpp"
#include "jgap/io/convert/QuipXmlConverter.hpp"
#include "jgap/utils/gap/StandardTabulation.hpp"

using namespace jgap;
using namespace jgap::validation;
namespace fs = std::filesystem;

static std::map<std::string, std::vector<double>> readTabGapH53b(const fs::path& h5_path) {
    std::map<std::string, std::vector<double>> result;
    HighFive::File file(h5_path.string(), HighFive::File::ReadOnly);
    for (const auto& name : file.listObjectNames()) {
        int dashes = 0;
        for (char c : name) if (c == '-') dashes++;
        if (dashes == 2) {
            HighFive::Group grp = file.getGroup(name);
            if (grp.exist("energies")) {
                std::vector<double> vals;
                grp.getDataSet("energies").read(vals);
                result[name] = vals;
            }
        }
    }
    return result;
}

static std::vector<double> readEamFloats(const fs::path& eam_path) {
    std::vector<double> vals;
    std::ifstream f(eam_path);
    if (!f) return vals;
    std::string line;
    // Skip 5 header lines
    for (int i = 0; i < 5 && std::getline(f, line); ++i);
    double v;
    while (f >> v) {
        vals.push_back(v);
    }
    return vals;
}

int main(int argc, char* argv[]) {
    ValidationReporter reporter("ValidateTabulation");

    fs::path res_dir = getValidationResourceDir();
    fs::path tables_dir = res_dir / "reference_tables";
    fs::path pots_dir = res_dir / "reference_pots";

    std::vector<std::string> run_names;
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--run" && i + 1 < argc) {
            run_names.push_back(argv[++i]);
        } else if (!arg.starts_with("-")) {
            run_names.push_back(arg);
        }
    }

    if (run_names.empty() && fs::exists(tables_dir)) {
        for (const auto& entry: fs::directory_iterator(tables_dir)) {
            if (entry.is_directory()) {
                std::string name = entry.path().filename().string();
                fs::path ref_pot = pots_dir / name / "gap.xml";
                fs::path ref_h5 = entry.path() / (name + ".tabgap.h5");
                if (fs::exists(ref_pot) && fs::exists(ref_h5)) {
                    run_names.push_back(name);
                }
            }
        }
        std::sort(run_names.begin(), run_names.end());
    }

    if (run_names.empty()) {
        std::cerr << "Error: No matching runs found in " << tables_dir << "\n";
        return 1;
    }

    std::cout << "======================================================================\n";
    std::cout << "JGAP Validation: Tabulation vs Reference tabGAP (Total Runs: " << run_names.size() << ")\n";
    std::cout << "  Tables Dir : " << tables_dir << "\n";
    std::cout << "  Pots Dir   : " << pots_dir << "\n";
    std::cout << "======================================================================\n";

    fs::path tmp_dir = fs::temp_directory_path() / "jgap_val_tabulation";
    fs::create_directories(tmp_dir);

    for (const auto& run_name: run_names) {
        fs::path ref_case_dir = tables_dir / run_name;
        fs::path ref_xml = pots_dir / run_name / "gap.xml";
        fs::path ref_h5 = ref_case_dir / (run_name + ".tabgap.h5");
        fs::path ref_eam = ref_case_dir / (run_name + ".eam.fs");

        std::cout << "\n>>> Tabulating and validating run: " << run_name << "\n";
        fs::path out_prefix = tmp_dir / ("tab_" + run_name);

        ValuePtr<Potential> potential = QuipXmlConverter::transform(ref_xml);
        if (!potential) {
            std::cerr << "Error: Failed to convert potential XML for " << run_name << "\n";
            reporter.recordFailure("Failed to convert potential XML for " + run_name);
            continue;
        }

        utils::StandardTabulationParams params;
        params.r_min_3b = 0.1;
        params.max_eam_density = 10.0;
        params.n_grid_2b = 1000;
        params.n_grid_3b = {20, 20, 20};

        utils::standardTabulation(*potential, out_prefix.string(), params);

        fs::path gen_h5 = tmp_dir / ("tab_" + run_name + ".tabgap.h5");
        fs::path gen_eam = tmp_dir / ("tab_" + run_name + ".eam.fs");

        if (!fs::exists(gen_h5)) {
            reporter.recordFailure("Generated .tabgap.h5 file missing for " + run_name);
            continue;
        }

        // 1. Compare 3B spline grid tables
        if (fs::exists(ref_h5)) {
            auto gen_3b = readTabGapH53b(gen_h5);
            auto ref_3b = readTabGapH53b(ref_h5);

            for (const auto& [triplet, q_vals] : ref_3b) {
                if (gen_3b.find(triplet) != gen_3b.end()) {
                    const auto& j_vals = gen_3b[triplet];
                    std::vector<double> j_mev, q_mev;
                    for (double v : j_vals) j_mev.push_back(v * 1e3);
                    for (double v : q_vals) q_mev.push_back(v * 1e3);

                    double rmse = computeRmse(j_mev, q_mev);
                    double nrmse = computeNrmsePct(j_mev, q_mev);
                    std::string test_id = run_name + "/3B_" + triplet;
                    reporter.assertLessOrEqual(test_id, "RMSE", rmse, 0.05, " meV");
                    reporter.assertLessOrEqual(test_id, "NRMSE", nrmse, 0.01, " %");
                }
            }
        }

        // 2. Compare EAM tables (.eam.fs)
        if (fs::exists(ref_eam) && fs::exists(gen_eam)) {
            auto gen_vals = readEamFloats(gen_eam);
            auto ref_vals = readEamFloats(ref_eam);
            if (!gen_vals.empty() && gen_vals.size() == ref_vals.size()) {
                std::vector<double> j_mev, q_mev;
                for (double v : gen_vals) j_mev.push_back(v * 1e3);
                for (double v : ref_vals) q_mev.push_back(v * 1e3);

                double rmse = computeRmse(j_mev, q_mev);
                double nrmse = computeNrmsePct(j_mev, q_mev);
                std::string test_id = run_name + "/EAM_fs";
                reporter.assertLessOrEqual(test_id, "RMSE", rmse, 0.01, " meV");
                reporter.assertLessOrEqual(test_id, "NRMSE", nrmse, 0.01, " %");
            }
        }
    }

    // Cleanup temp files
    try {
        fs::remove_all(tmp_dir);
    } catch (...) {}

    return reporter.summarizeAndExit();
}
