#include <cmath>
#include <filesystem>
#include <iostream>
#include <string>
#include <vector>

#include "common/ValidationUtils.hpp"

#include "jgap/core/UnseqFor.hpp"
#include "jgap/core/atomic/Atoms.hpp"
#include "jgap/core/potentials/tabgap/TabGapPotential.hpp"
#include "jgap/io/tabgap/TabGapIO.hpp"

using namespace jgap;
using namespace jgap::validation;
namespace fs = std::filesystem;

int main(int argc, char* argv[]) {
    ValidationReporter reporter("ValidateTabGapEnergyEval");

    fs::path res_dir = getValidationResourceDir();
    fs::path root = findRepoRoot();

    fs::path test_xyz = root / "test" / "resources" / "structure-databases" / "feni-test.xyz";
    fs::path lammps_dir = res_dir / "reference_lammps";
    fs::path tables_dir = res_dir / "reference_tables";

    std::vector<std::string> run_names;
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--run" && i + 1 < argc) {
            run_names.push_back(argv[++i]);
        } else if (!arg.starts_with("-")) {
            run_names.push_back(arg);
        }
    }

    if (run_names.empty() && fs::exists(lammps_dir)) {
        for (const auto& entry: fs::directory_iterator(lammps_dir)) {
            if (entry.is_directory()) {
                std::string name = entry.path().filename().string();
                fs::path h5 = tables_dir / name / (name + ".tabgap.h5");
                if (fs::exists(entry.path() / "pred.xyz") && fs::exists(h5)) {
                    run_names.push_back(name);
                }
            }
        }
        std::sort(run_names.begin(), run_names.end());
    }

    if (run_names.empty()) {
        std::cerr << "Error: No matching runs found in " << lammps_dir << "\n";
        return 1;
    }

    std::cout << "======================================================================\n";
    std::cout << "JGAP Validation: TabGapPotential vs LAMMPS (Total Runs: " << run_names.size() << ")\n";
    std::cout << "  Test XYZ     : " << test_xyz << "\n";
    std::cout << "  LAMMPS Dir   : " << lammps_dir << "\n";
    std::cout << "======================================================================\n";

    for (const auto& run_name: run_names) {
        fs::path table_dir = tables_dir / run_name;
        fs::path h5_path = table_dir / (run_name + ".tabgap.h5");
        fs::path eam_path = table_dir / (run_name + ".eam.fs");
        fs::path ref_lammps_xyz = lammps_dir / run_name / "pred.xyz";

        std::cout << "\n>>> Validating run: " << run_name << "\n";

        std::vector<std::string> pot_files = {h5_path.string()};
        if (fs::exists(eam_path)) {
            pot_files.push_back(eam_path.string());
        }

        TabGapPotential potential = TabGapIO::read(pot_files);

        std::vector<Atoms> frames = Atoms::readAtoms(test_xyz.string());
        unseqForEach(frames.begin(), frames.end(), [&](Atoms& atoms) {
            atoms.setEnergyAndDerivatives(potential.calculateEnergy(atoms));
        });

        MainXYZPropertyNames prop_names;
        prop_names.virials = "virial";
        const std::vector<Atoms> ref_frames = Atoms::readAtoms(ref_lammps_xyz.string(), prop_names);

        if (frames.size() != ref_frames.size()) {
            std::cerr << "Error: Frame count mismatch (" << frames.size() << " vs " << ref_frames.size() << ") for run " << run_name << "\n";
            reporter.recordFailure("Frame count mismatch for " + run_name);
            continue;
        }

        std::vector<double> jgap_energies_mev_per_atom, ref_energies_mev_per_atom;
        std::vector<double> jgap_forces_mev, ref_forces_mev;

        for (size_t f = 0; f < frames.size(); ++f) {
            const auto& j_atoms = frames[f];
            const auto& q_atoms = ref_frames[f];
            size_t n = j_atoms.nAtoms();
            if (n == 0) continue;

            double e_j = j_atoms.getEnergy().value_or(0.0) / static_cast<double>(n) * 1e3;
            double e_q = q_atoms.getEnergy().value_or(0.0) / static_cast<double>(n) * 1e3;
            jgap_energies_mev_per_atom.push_back(e_j);
            ref_energies_mev_per_atom.push_back(e_q);

            const auto fj_opt = j_atoms.getForces();
            const auto fq_opt = q_atoms.getForces();
            if (fj_opt.has_value() && fq_opt.has_value()) {
                const auto& fj = fj_opt.value();
                const auto& fq = fq_opt.value();
                for (size_t i = 0; i < fj.size(); ++i) {
                    jgap_forces_mev.push_back(fj[i].x * 1e3);
                    jgap_forces_mev.push_back(fj[i].y * 1e3);
                    jgap_forces_mev.push_back(fj[i].z * 1e3);
                    ref_forces_mev.push_back(fq[i].x * 1e3);
                    ref_forces_mev.push_back(fq[i].y * 1e3);
                    ref_forces_mev.push_back(fq[i].z * 1e3);
                }
            }
        }

        double e_rmse = computeRmse(jgap_energies_mev_per_atom, ref_energies_mev_per_atom);
        double f_rmse = computeRmse(jgap_forces_mev, ref_forces_mev);

        std::cout << "--- " << run_name << " TabGap vs LAMMPS Metrics ---\n";
        reporter.assertLessOrEqual(run_name + "/Energy", "RMSE", e_rmse, 0.01, " meV/atom");
        reporter.assertLessOrEqual(run_name + "/Forces", "RMSE", f_rmse, 0.50, " meV/Å");
    }

    return reporter.summarizeAndExit();
}
