#include <cmath>
#include <filesystem>
#include <iostream>
#include <string>
#include <vector>

#include "common/ValidationUtils.hpp"

#include "jgap/core/UnseqFor.hpp"
#include "jgap/core/atomic/Atoms.hpp"
#include "jgap/io/convert/QuipXmlConverter.hpp"
#include "jgap/serialization/SerializationRegistry.hpp"

using namespace jgap;
using namespace jgap::validation;
namespace fs = std::filesystem;

int main(int argc, char* argv[]) {
    ValidationReporter reporter("ValidateGapEnergyEval");

    fs::path res_dir = getValidationResourceDir();
    fs::path root = findRepoRoot();

    fs::path test_xyz = root / "test" / "resources" / "structure-databases" / "feni-test.xyz";
    fs::path preds_dir = res_dir / "reference_preds";
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

    if (run_names.empty() && fs::exists(preds_dir)) {
        for (const auto& entry: fs::directory_iterator(preds_dir)) {
            if (entry.is_directory()) {
                std::string name = entry.path().filename().string();
                if (fs::exists(entry.path() / "pred.xyz") && fs::exists(pots_dir / name / "gap.xml")) {
                    run_names.push_back(name);
                }
            }
        }
        std::sort(run_names.begin(), run_names.end());
    }

    if (run_names.empty()) {
        std::cerr << "Error: No matching runs found in " << preds_dir << "\n";
        return 1;
    }

    std::cout << "======================================================================\n";
    std::cout << "JGAP Validation: Test Set Predictions vs QUIP (Total Runs: " << run_names.size() << ")\n";
    std::cout << "  Test XYZ     : " << test_xyz << "\n";
    std::cout << "  Reference Dir: " << preds_dir << "\n";
    std::cout << "======================================================================\n";

    for (const auto& run_name: run_names) {
        fs::path pot_xml = pots_dir / run_name / "gap.xml";
        fs::path ref_pred_xyz = preds_dir / run_name / "pred.xyz";

        std::cout << "\n>>> Validating run: " << run_name << "\n";

        ValuePtr<Potential> potential = QuipXmlConverter::transform(pot_xml);
        if (!potential) {
            std::cerr << "Error: Failed to convert potential XML for " << run_name << "\n";
            reporter.recordFailure("Failed to convert potential XML for " + run_name);
            continue;
        }

        std::vector<Atoms> frames = Atoms::readAtoms(test_xyz.string());
        unseqForEach(frames.begin(), frames.end(), [&](Atoms& atoms) {
            atoms.setEnergyAndDerivatives(potential->calculateEnergy(atoms));
        });

        MainXYZPropertyNames prop_names;
        prop_names.virials = "virial";
        const std::vector<Atoms> ref_frames = Atoms::readAtoms(ref_pred_xyz.string(), prop_names);

        if (frames.size() != ref_frames.size()) {
            std::cerr << "Error: Mismatch in number of frames (" << frames.size()
                      << " vs " << ref_frames.size() << ") for run " << run_name << "\n";
            reporter.recordFailure("Mismatch in number of frames for " + run_name);
            continue;
        }

    std::vector<double> jgap_energies_per_atom, ref_energies_per_atom;
    std::vector<double> jgap_forces, ref_forces;
    std::vector<double> jgap_virials_per_atom, ref_virials_per_atom;

    for (size_t f = 0; f < frames.size(); ++f) {
        const auto& j_atoms = frames[f];
        const auto& q_atoms = ref_frames[f];
        size_t n = j_atoms.nAtoms();
        if (n == 0) continue;

        // Energy per atom (eV -> meV)
        double e_j = j_atoms.getEnergy().value_or(0.0) / static_cast<double>(n) * 1e3;
        double e_q = q_atoms.getEnergy().value_or(0.0) / static_cast<double>(n) * 1e3;
        jgap_energies_per_atom.push_back(e_j);
        ref_energies_per_atom.push_back(e_q);

        // Forces (eV/Å -> meV/Å)
        const auto fj_opt = j_atoms.getForces();
        const auto fq_opt = q_atoms.getForces();
        if (fj_opt.has_value() && fq_opt.has_value()) {
            const auto& fj = fj_opt.value();
            const auto& fq = fq_opt.value();
            for (size_t i = 0; i < fj.size(); ++i) {
                jgap_forces.push_back(fj[i].x * 1e3);
                jgap_forces.push_back(fj[i].y * 1e3);
                jgap_forces.push_back(fj[i].z * 1e3);
                ref_forces.push_back(fq[i].x * 1e3);
                ref_forces.push_back(fq[i].y * 1e3);
                ref_forces.push_back(fq[i].z * 1e3);
            }
        }

        // Virials (eV -> meV per atom)
        const auto vj_opt = j_atoms.getVirials();
        const auto vq_opt = q_atoms.getVirials();
        if (vj_opt.has_value() && vq_opt.has_value()) {
            const auto& vj = vj_opt.value();
            const auto& vq = vq_opt.value();
            double inv_n = 1e3 / static_cast<double>(n);
            jgap_virials_per_atom.push_back(vj.xx * inv_n);
            jgap_virials_per_atom.push_back(vj.xy * inv_n);
            jgap_virials_per_atom.push_back(vj.xz * inv_n);
            jgap_virials_per_atom.push_back(vj.yy * inv_n);
            jgap_virials_per_atom.push_back(vj.yz * inv_n);
            jgap_virials_per_atom.push_back(vj.zz * inv_n);

            ref_virials_per_atom.push_back(vq.xx * inv_n);
            ref_virials_per_atom.push_back(vq.xy * inv_n);
            ref_virials_per_atom.push_back(vq.xz * inv_n);
            ref_virials_per_atom.push_back(vq.yy * inv_n);
            ref_virials_per_atom.push_back(vq.yz * inv_n);
            ref_virials_per_atom.push_back(vq.zz * inv_n);
        }
    }

    double e_rmse = computeRmse(jgap_energies_per_atom, ref_energies_per_atom);
    double f_rmse = computeRmse(jgap_forces, ref_forces);
    double v_rmse = computeRmse(jgap_virials_per_atom, ref_virials_per_atom);

    double e_nrmse = computeNrmsePct(jgap_energies_per_atom, ref_energies_per_atom);
    double f_nrmse = computeNrmsePct(jgap_forces, ref_forces);
    double v_nrmse = computeNrmsePct(jgap_virials_per_atom, ref_virials_per_atom);

        std::cout << "--- " << run_name << " Metrics vs QUIP ---\n";
        reporter.assertLessOrEqual(run_name + "/Energy", "RMSE", e_rmse, 0.01, " meV/atom");
        reporter.assertLessOrEqual(run_name + "/Forces", "RMSE", f_rmse, 0.01, " meV/Å");
        reporter.assertLessOrEqual(run_name + "/Virials", "RMSE", v_rmse, 0.01, " meV/atom");

        reporter.assertLessOrEqual(run_name + "/Energy", "NRMSE", e_nrmse, 0.001, " %");
        reporter.assertLessOrEqual(run_name + "/Forces", "NRMSE", f_nrmse, 0.001, " %");
        reporter.assertLessOrEqual(run_name + "/Virials", "NRMSE", v_nrmse, 0.001, " %");
    }

    return reporter.summarizeAndExit();
}
