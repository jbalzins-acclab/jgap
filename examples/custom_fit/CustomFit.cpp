// Example: a super-custom multi-element Fe-Ni GAP demonstrating advanced & experimental components:
//   * Multiple kernel types: WendlandKernel (2b), CauchyKernel (EAM), SquaredExpKernel (MEAM);
//   * MEAM 3-body descriptor: ThreeBodySum with MeamTransformation (Legendre polynomial expansion);
//   * Polynomial Cutoff: PerriotPolynomialCutoff (smooth C2 polynomial cutoff function);
//   * Hybrid EAM density functions: FSGen (Fe-Fe), Coscutoff (Ni-Ni), Polycutoff (Fe-Ni);
//   * Adaptive Regularization: ScaledRegularizationRules (scales sigmas with max atomic force);
//   * Fit: ElementIncrementalQRGapFit (streaming out-of-core QR with configurable RAM limit).
//
// Usage: custom_fit [training.xyz] [output_prefix] [--ram-limit <gb>]
//   defaults: test/structure-databases/feni-train.xyz  feni-custom

#include <chrono>
#include <filesystem>
#include <iostream>
#include <set>
#include <string>
#include <vector>

#include "jgap/core/UnseqFor.hpp"
#include "jgap/core/fit/gap/regularization/PerConfigTypeRegularizationRules.hpp"
#include "jgap/core/fit/gap/regularization/ScaledRegularizationRules.hpp"
#include "jgap/core/transform/manybody/TwoBodySum.hpp"
#include "jgap/impl/cutoff/PerriotPolynomialCutoff.hpp"
#include "jgap/impl/fit/gap/ElementIncrementalQRGapFit.hpp"
#include "jgap/impl/kernels/CauchyKernel.hpp"
#include "jgap/impl/kernels/SquaredExpKernel.hpp"
#include "jgap/impl/kernels/WendlandKernel.hpp"
#include "jgap/impl/transform/manybody/ThreeBodySum.hpp"
#include "jgap/impl/transform/nbody/3b/MeamTransformation.hpp"
#include "jgap/jgap.hpp"

using namespace jgap;
using namespace jgap::utils;

namespace {
    // ---- baseline parameters ----
    constexpr size_t SEED = 120;

    constexpr double EAM_CUTOFF = 5.0;
    constexpr double EAM_ENERGY_SCALE = 1.0;
    constexpr size_t EAM_N_SPARSE = 20;

    constexpr double MEAM_CUTOFF = 4.0;
    constexpr double MEAM_WIDTH = 0.6;
    constexpr double MEAM_ENERGY_SCALE = 1.0;
    constexpr size_t MEAM_N_SPARSE = 50;

    constexpr double CUTOFF_2B = 5.0;
    constexpr double WIDTH_2B = 0.5;
    constexpr double ENERGY_SCALE_2B = 10.0;
    constexpr size_t N_SPARSE_2B = 20;

    bool isNi(const Species& s) { return s.symbol() == "Ni"; }

    /// Picks a different EAM pair (density) function per element pair, demonstrating hybrid density types.
    ValuePtr<EamPairFunction> makeEamPairFunction(const Species& central, const Species& contributor) {
        const int n_ni = (isNi(central) ? 1 : 0) + (isNi(contributor) ? 1 : 0);
        if (n_ni == 0) {
            return FSGenPairFunction(EAM_CUTOFF, /*degree=*/3.0); // Fe-Fe
        }
        if (n_ni == 2) {
            return CoscutoffPairFunction(EAM_CUTOFF, /*r_min=*/0.0); // Ni-Ni
        }
        return PolycutoffPairFunction(EAM_CUTOFF, /*r_min=*/0.0); // mixed Fe-Ni
    }
}

int main(int argc, char** argv) {
    CurrentLogger::initDefault({.stdout_log_debug = true});

    std::vector<std::string> positional_args;
    double ram_limit_gb = 2.0;

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--ram-limit" && i + 1 < argc) {
            ram_limit_gb = std::stod(argv[++i]);
        } else if (arg.starts_with("--ram-limit=")) {
            ram_limit_gb = std::stod(arg.substr(12));
        } else if (arg.starts_with("ram_limit=")) {
            ram_limit_gb = std::stod(arg.substr(10));
        } else if (arg == "-h" || arg == "--help") {
            std::cout << "Usage: " << argv[0] << " [training.xyz] [output_prefix] [--ram-limit <gb>]\n"
                      << "  --ram-limit <gb>   RAM limit in GB for ElementIncrementalQRGapFit (default: 2.0)\n";
            return 0;
        } else {
            positional_args.push_back(arg);
        }
    }

    std::string default_training = "test/structure-databases/feni-train.xyz";
    if (!std::filesystem::exists(default_training) && std::filesystem::exists("../../test/structure-databases/feni-train.xyz")) {
        default_training = "../../test/structure-databases/feni-train.xyz";
    }

    const std::string training_file = !positional_args.empty() ? positional_args[0] : default_training;
    const std::string output_prefix = positional_args.size() > 1 ? positional_args[1] : "feni-custom";

    const auto start = std::chrono::steady_clock::now();

    JGAP_LOG_INFO("Experimental Custom FeNi fit on {} using ElementIncrementalQRGapFit (RAM limit: {} GB)", training_file, ram_limit_gb);
    const auto training_data = Atoms::readAtoms(training_file);

    // Elements present in the training data.
    std::set<Species> elements;
    for (const auto& atoms: training_data) {
        for (const auto& s: atoms.getSpecies()) {
            elements.insert(s);
        }
    }

    GapPotential potential;

    // ===== 1. EAM: CauchyKernel (rational quadratic) + Hybrid pair functions =====
    const auto eam_kernel = CauchyKernel<1, 0>(EAM_ENERGY_SCALE, {1.0});
    const HistogramUniformSparsifier<1> eam_sparsifier(SEED, EAM_N_SPARSE);
    for (const Species& central: elements) {
        auto aggregator = std::make_unique<TwoBodySum<1>>(central);
        for (const Species& contributor: elements) {
            aggregator->extend({central, contributor}, makeEamPairFunction(central, contributor));
        }
        potential.addComponent(ManyBodyGapComponent(
            ValuePtr<NBodyAggregator<1>>(std::move(aggregator)), eam_kernel, eam_sparsifier, training_data
        ));
    }

    // ===== 2. MEAM: ThreeBodySum + MeamTransformation (Legendre expansion) + SquaredExpKernel =====
    const ValuePtr<CutoffFunction> meam_cutoff = PerriotPolynomialCutoff(MEAM_CUTOFF, MEAM_WIDTH);
    const ValuePtr<MeamTransformation> meam_trans = MeamTransformation(meam_cutoff);
    const auto meam_kernel = SquaredExpKernel<3, 0>(MEAM_ENERGY_SCALE, {1.0, 1.0, 1.0});
    const HistogramUniformSparsifier<3> meam_sparsifier(SEED, MEAM_N_SPARSE);

    for (const Species& central: elements) {
        auto meam_aggregator = std::make_unique<ThreeBodySum<3>>(central);
        for (const Species& c1: elements) {
            for (const Species& c2: elements) {
                meam_aggregator->extend({central, c1, c2}, meam_trans);
            }
        }
        potential.addComponent(ManyBodyGapComponent<3, SquaredExpKernel<3, 0>>(
            ValuePtr<NBodyAggregator<3>>(std::move(meam_aggregator)), meam_kernel, meam_sparsifier, training_data
        ));
    }

    // ===== 3. 2-body: WendlandKernel (compact support) + PerriotPolynomialCutoff =====
    const ValuePtr<TwoBodyTransformation<2>> trans2 =
        PairDistanceTransformation(PerriotPolynomialCutoff(CUTOFF_2B, WIDTH_2B));
    const auto kernel2 = WendlandKernel<1, 1>(ENERGY_SCALE_2B, {1.0});
    const HistogramUniformSparsifier<2> sparsifier2(SEED, N_SPARSE_2B, std::array{true, false});
    potential.addComponents(
        createTwoBodyComponents<2, WendlandKernel<1, 1>>(training_data, trans2, kernel2, sparsifier2)
    );

    // ===== 4. External potentials: isolated-atom energies + ScreenedCoulomb repulsion =====
    potential.optional_external_potential = CompositePotential{{
        {"isolated", IsolatedAtomPotential{training_data}},
        {"screened_coulomb", ScreenedCoulombPotential{training_data}},
    }};

    // ===== 5. Regularization: ScaledRegularizationRules (force-adaptive scaling) =====
    const PerConfigTypeRegularizationRules base_rules(
        PerConfigTypeSigmas(0.002, 0.1, 0.2),
        "isolated_atom:0.0001:0.04:0.04:0.0:liquid:0.01:0.5:2.0:0.0:dimer:0.01:0.5:2.0:0.0:"
        "short_range:0.01:0.5:2.0:0.0:liquid_surface_100:0.01:0.5:2.0:0.0:"
        "liquid_surface_110:0.01:0.5:2.0:0.0:liquid_surface_111:0.01:0.5:2.0:0.0:"
        "gamma_surface:0.002:0.08:0.5:0.0:liquid_high:0.02:0.8:5.0:0.0:"
        "binary_alloy_melting:0.01:0.5:2.0:0.0:binary_alloy_short_range:0.01:0.5:2.0:0.0"
    );
    const ScaledRegularizationRules regularization(base_rules, /*force_scale=*/0.5, /*min_scale=*/1.0);

    // ===== 6. Out-of-core streaming QR fit =====
    ElementIncrementalQRGapFit fitter(1e-8, ram_limit_gb);
    auto sigmas = regularization.determineForAll(training_data);
    fitter.fit(potential, training_data, sigmas);

    const std::string potential_file = output_prefix + ".jgap.h5";
    SerializationRegistry<Potential>::serialize(potential, potential_file);
    JGAP_LOG_INFO("Saved fitted experimental potential to {}", potential_file);

    // Verify evaluation directly on training structure
    if (!training_data.empty()) {
        auto eval_result = potential.calculateEnergy(training_data[0]);
        JGAP_LOG_INFO("Sample prediction on frame 0: energy = {:.6f} eV (ref = {:.6f} eV)",
                      eval_result.value, training_data[0].getEnergy().value_or(0.0));
    }

    std::cout << "Execution time: " << formatDuration(elapsedMillisSince(start)) << std::endl;
    return 0;
}
