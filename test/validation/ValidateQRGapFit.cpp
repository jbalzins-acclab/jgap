#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <numeric>
#include <string>
#include <vector>

#include "common/ValidationUtils.hpp"

#include "jgap/core/atomic/Atoms.hpp"
#include "jgap/core/fit/gap/regularization/PerConfigTypeRegularizationRules.hpp"
#include "jgap/core/potentials/gap/GapPotential.hpp"
#include "jgap/impl/fit/gap/QRGapFit.hpp"
#include "jgap/serialization/SerializationRegistry.hpp"

#include "jgap/io/convert/QuipXmlConverter.hpp"

using namespace jgap;
using namespace jgap::validation;
namespace fs = std::filesystem;

static std::vector<double> getAllCoefficients(const GapPotential& pot) {
    std::vector<double> coeffs;
    for (const auto& comp : pot.getComponents()) {
        const auto& c = comp->getCoefficients();
        coeffs.insert(coeffs.end(), c.begin(), c.end());
    }
    return coeffs;
}

int main(int argc, char* argv[]) {
    ValidationReporter reporter("ValidateQRGapFit");

    fs::path res_dir = getValidationResourceDir();
    fs::path pots_dir = res_dir / "reference_pots";

    std::vector<std::string> run_names;
    bool quick_mode = false;
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--run" && i + 1 < argc) {
            run_names.push_back(argv[++i]);
        } else if (arg == "--quick") {
            quick_mode = true;
        } else if (!arg.starts_with("-")) {
            run_names.push_back(arg);
        }
    }

    if (quick_mode && run_names.empty()) {
        run_names.push_back("feni_200_100_0");
    }

    if (run_names.empty() && fs::exists(pots_dir)) {
        for (const auto& entry: fs::directory_iterator(pots_dir)) {
            if (entry.is_directory()) {
                std::string name = entry.path().filename().string();
                if (fs::exists(entry.path() / "gap.xml") && fs::exists(entry.path() / "train.xyz")) {
                    run_names.push_back(name);
                }
            }
        }
        std::sort(run_names.begin(), run_names.end());
    }

    if (run_names.empty()) {
        std::cerr << "Error: No matching reference potential runs found in " << pots_dir << "\n";
        return 1;
    }

    std::cout << "======================================================================\n";
    std::cout << "JGAP Validation: Fit Coefficients vs Reference QUIP Potentials (Total Runs: " << run_names.size() << ")\n";
    std::cout << "  Reference Dir: " << pots_dir << "\n";
    std::cout << "======================================================================\n";

    const PerConfigTypeRegularizationRules regularization(
        PerConfigTypeSigmas(0.002, 0.1, 0.2),
        "isolated_atom:0.0001:0.04:0.04:0.0:"
        "liquid:0.01:0.5:2.0:0.0:"
        "dimer:0.01:0.5:2.0:0.0:"
        "short_range:0.01:0.5:2.0:0.0:"
        "liquid_surface_100:0.01:0.5:2.0:0.0:"
        "liquid_surface_110:0.01:0.5:2.0:0.0:"
        "liquid_surface_111:0.01:0.5:2.0:0.0:"
        "gamma_surface:0.002:0.08:0.5:0.0:"
        "liquid_high:0.02:0.8:5.0:0.0:"
        "binary_alloy_melting:0.01:0.5:2.0:0.0:"
        "binary_alloy_short_range:0.01:0.5:2.0:0.0"
    );

    for (const auto& run_name: run_names) {
        fs::path ref_case_dir = pots_dir / run_name;
        fs::path ref_xml = ref_case_dir / "gap.xml";
        fs::path train_xyz = ref_case_dir / "train.xyz";

        std::cout << "\n>>> Validating run: " << run_name << "\n";

        ValuePtr<Potential> ref_pot_ptr = QuipXmlConverter::transform(ref_xml);
        auto* ref_gap = dynamic_cast<GapPotential*>(ref_pot_ptr.get());
        if (!ref_gap) {
            std::cerr << "Error: Failed to convert reference potential as GapPotential for " << run_name << ".\n";
            reporter.recordFailure("Failed to convert reference potential for " + run_name);
            continue;
        }

        MainXYZPropertyNames prop_names;
        prop_names.virials = "virial_fit";
        const std::vector<Atoms> train_data = Atoms::readAtoms(train_xyz.string(), prop_names);
        auto sigmas = regularization.determineForAll(train_data);

        std::cout << "Fitting JGAP potential using reference sparse points..." << std::flush;
        GapPotential fit_pot(*ref_gap);
        QRGapFit fitter(1e-8);
        fitter.fit(fit_pot, train_data, sigmas);
        std::cout << " Done.\n";

        auto coeffs_ref = getAllCoefficients(*ref_gap);
        auto coeffs_jgap = getAllCoefficients(fit_pot);

        double overall_nrmse = computeNrmsePct(coeffs_jgap, coeffs_ref);
        double overall_sigrel = computeSigRelPct(coeffs_jgap, coeffs_ref);
        double overall_cosine = computeCosineSim(coeffs_jgap, coeffs_ref);

        reporter.assertLessOrEqual(run_name + "/Coefficients_Overall", "NRMSE", overall_nrmse, 0.01, "%");
        reporter.assertLessOrEqual(run_name + "/Coefficients_Overall", "SigRel", overall_sigrel, 0.05, "%");
        reporter.assertGreaterOrEqual(run_name + "/Coefficients_Overall", "CosineSim", overall_cosine, 0.999999);

        for (size_t c = 0; c < fit_pot.getComponents().size(); ++c) {
            auto q_c = ref_gap->getComponents()[c]->getCoefficients();
            auto j_c = fit_pot.getComponents()[c]->getCoefficients();
            double c_nrmse = computeNrmsePct(j_c, q_c);
            double c_cosine = computeCosineSim(j_c, q_c);
            std::string comp_name = run_name + "/Component_" + std::to_string(c);
            reporter.assertLessOrEqual(comp_name, "NRMSE", c_nrmse, 0.01, "%");
            reporter.assertGreaterOrEqual(comp_name, "CosineSim", c_cosine, 0.999999);
        }
    }

    return reporter.summarizeAndExit();
}
