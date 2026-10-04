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
#include "jgap/ext/fit/gap/QRGapFit.hpp"
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
    fs::path ref_case_dir = res_dir / "reference_pots" / "feni_200_100_0";

    if (!fs::exists(ref_case_dir)) {
        std::cerr << "Error: Reference potential directory not found: " << ref_case_dir << "\n";
        return 1;
    }

    fs::path ref_xml = ref_case_dir / "gap.xml";
    fs::path train_xyz = ref_case_dir / "train.xyz";
    if (!fs::exists(train_xyz)) {
        train_xyz = findRepoRoot() / "test" / "resources" / "structure-databases" / "feni-train.xyz";
    }

    std::cout << "======================================================================\n";
    std::cout << "JGAP Validation: Fit Coefficients vs Reference QUIP Potential\n";
    std::cout << "  Reference Potential : " << ref_xml << "\n";
    std::cout << "  Training Set        : " << train_xyz << "\n";
    std::cout << "======================================================================\n";

    ValuePtr<Potential> ref_pot_ptr = QuipXmlConverter::transform(ref_xml);
    auto* ref_gap = dynamic_cast<GapPotential*>(ref_pot_ptr.get());
    if (!ref_gap) {
        std::cerr << "Error: Failed to convert reference potential as GapPotential.\n";
        return 1;
    }

    MainXYZPropertyNames prop_names;
    prop_names.virials = "virial_fit";
    const std::vector<Atoms> train_data = Atoms::readAtoms(train_xyz.string(), prop_names);

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

    std::cout << "\n--- Coefficient Comparison Results ---\n";
    reporter.assertLessOrEqual("Coefficients_Overall", "NRMSE", overall_nrmse, 0.01, "%");
    reporter.assertLessOrEqual("Coefficients_Overall", "SigRel", overall_sigrel, 0.05, "%");
    reporter.assertGreaterOrEqual("Coefficients_Overall", "CosineSim", overall_cosine, 0.999999);

    for (size_t c = 0; c < fit_pot.getComponents().size(); ++c) {
        auto q_c = ref_gap->getComponents()[c]->getCoefficients();
        auto j_c = fit_pot.getComponents()[c]->getCoefficients();
        double c_nrmse = computeNrmsePct(j_c, q_c);
        double c_cosine = computeCosineSim(j_c, q_c);
        std::string comp_name = "Component_" + std::to_string(c);
        reporter.assertLessOrEqual(comp_name, "NRMSE", c_nrmse, 0.01, "%");
        reporter.assertGreaterOrEqual(comp_name, "CosineSim", c_cosine, 0.999999);
    }

    return reporter.summarizeAndExit();
}
