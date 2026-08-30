// Example: a more flexible FeNi fit built by hand (instead of standardGapFit), showing how to vary the
// per-component setup:
//   * EAM: a different pair (density) function type per element pair — FSGen for Fe-Fe, Coscutoff for
//     Ni-Ni, Polycutoff for the mixed pairs;
//   * 3-body: a per-triplet energy scale that grows slightly with the number of Ni atoms in the triplet;
//   * 2-body / external (isolated + ZBL): as in the standard fit.
//
// Baseline numbers (cutoffs, sparse counts, scales, regularization) are taken from the old parameter
// file test-local/fit-params.json. jgap has no JSON reader, so they are inlined here as clearly-labelled
// constants rather than parsed.
//
// Usage: custom_fit [training.xyz] [output_prefix]
//   defaults: resources/xyz-samples/feni-train.xyz  feni-custom

#include <chrono>
#include <iostream>
#include <jgap/ext/fit/gap/QRGapFit.hpp>
#include <set>
#include <string>

#include "jgap/core/UnseqFor.hpp"
#include "jgap/core/atomic/species/composition/Species3AtomicSorted.hpp"
#include "jgap/core/fit/gap/regularization/PerConfigTypeRegularizationRules.hpp"
#include "jgap/experimental/fit/gap/StreamingQrGapFit.hpp"
#include "jgap/jgap.hpp"

using namespace jgap;
using namespace jgap::utils;

namespace {
    // ---- baseline parameters, from test-local/fit-params.json ----
    constexpr size_t SEED = 120;

    constexpr Real EAM_CUTOFF = 5.0; // json eam pair_function cutoff
    constexpr Real EAM_ENERGY_SCALE = 1.0; // json eam kernel_setup energy_scale
    constexpr size_t EAM_N_SPARSE = 20; // json eam kernel_setup n_sparse

    constexpr Real CUTOFF_3B = 4.0; // json 3b cutoff
    constexpr Real WIDTH_3B = 0.6; // json 3b cutoff_transition_width
    constexpr Real ENERGY_SCALE_3B = 1.0; // json 3b energy_scale (base, before the Ni adjustment)
    constexpr size_t N_SPARSE_3B = 500; // json 3b n_sparse (per triplet — reduce for a quick run)

    constexpr Real CUTOFF_2B = 5.0; // json 2b cutoff
    constexpr Real WIDTH_2B = 0.5; // json 2b cutoff_transition_width
    constexpr Real ENERGY_SCALE_2B = 10.0; // json 2b energy_scale
    constexpr size_t N_SPARSE_2B = 20; // json 2b n_sparse

    bool isNi(const Species& s) { return s.symbol() == "Ni"; }

    /// Picks a different EAM pair (density) function per element pair, to exercise all three kinds.
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

    const std::string training_file = argc > 1 ? argv[1] : "../../resources/xyz-samples/feni-train.xyz";
    const std::string output_prefix = argc > 2 ? argv[2] : "feni-lsmr";

    const auto start = std::chrono::steady_clock::now();

    JGAP_LOG_INFO("Custom FeNi fit on {}", training_file);
    const auto training_data = Atoms::readAtoms(training_file);

    // Elements present in the training data.
    std::set<Species> elements;
    for (const auto& atoms: training_data) {
        for (const auto& s: atoms.getSpecies()) {
            elements.insert(s);
        }
    }

    GapPotential potential;

    // ===== EAM: one component per central species, with a per-pair density function =====
    const auto eam_kernel = SquaredExpKernel<1, 0>(EAM_ENERGY_SCALE, {1.0});
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

    // ===== 3-body: per-triplet energy scale (grows slightly with the number of Ni atoms in the triplet) =====
    const ValuePtr<ThreeBodyTransformation<4>> trans3 = Angle3bTransformation(CosCutoff(CUTOFF_3B, WIDTH_3B));
    const HistogramUniformSparsifier<4> sparsifier3(SEED, N_SPARSE_3B, std::array{true, true, true, false});

    std::set<Species3AtomicSorted> triplets;
    for (const auto& atoms: training_data) {
        NeighbourLists nl(atoms, trans3->getCutoffs().maxOverall());
        auto sets = Species3AtomicSorted::getAll(nl);
        triplets.insert(sets.begin(), sets.end());
    }
    for (const auto& triplet: triplets) {
        const int n_ni =
            (isNi(triplet.root) ? 1 : 0) + (isNi(triplet.nodes[0]) ? 1 : 0) + (isNi(triplet.nodes[1]) ? 1 : 0);
        const Real energy_scale = ENERGY_SCALE_3B * (1.0 + 0.1 * static_cast<Real>(n_ni)); // mild Ni boost
        const auto kernel3 = SquaredExpKernel<3, 1>(energy_scale, {1.0, 1.0, 1.0});

        potential.addComponent(
            ThreeBodyGapComponent<4, SquaredExpKernel<3, 1>>(triplet, trans3, kernel3, sparsifier3, training_data)
        );
    }

    // ===== 2-body: one shared setup (as in the standard fit) =====
    const ValuePtr<TwoBodyTransformation<2>> trans2 = PairDistanceTransformation(CosCutoff(CUTOFF_2B, WIDTH_2B));
    const auto kernel2 = SquaredExpKernel<1, 1>(ENERGY_SCALE_2B, {1.0});
    const HistogramUniformSparsifier<2> sparsifier2(SEED, N_SPARSE_2B, std::array{true, false});
    potential.addComponents(
        TwoBodyGapComponent<2, SquaredExpKernel<1, 1>>::createComponents(training_data, trans2, kernel2, sparsifier2)
    );

    // ===== external: isolated-atom energies + ScreenedCoulomb repulsion =====
    potential.optional_external_potential = CompositePotential{{
        {"isolated", IsolatedAtomPotential{training_data}},
        {"screened_coulomb", ScreenedCoulombPotential{training_data}},
    }};

    // ===== setup regularization =====
    const PerConfigTypeRegularizationRules regularization(
        PerConfigTypeSigmas(0.002, 0.1, 0.2),
        "isolated_atom:0.0001:0.04:0.04:0.0:liquid:0.01:0.5:2.0:0.0:dimer:0.01:0.5:2.0:0.0:"
        "short_range:0.01:0.5:2.0:0.0:liquid_surface_100:0.01:0.5:2.0:0.0:"
        "liquid_surface_110:0.01:0.5:2.0:0.0:liquid_surface_111:0.01:0.5:2.0:0.0:"
        "gamma_surface:0.002:0.08:0.5:0.0:liquid_high:0.02:0.8:5.0:0.0:"
        "binary_alloy_melting:0.01:0.5:2.0:0.0:binary_alloy_short_range:0.01:0.5:2.0:0.0"
    );

    // ===== fit =====
    StreamingQrGapFit fitter(1e-8, 1000);
    auto sigmas = regularization.determineForAll(training_data);
    fitter.fit(potential, training_data, sigmas);

    const std::string potential_file = output_prefix + ".jgap.h5";
    SerializationRegistry<Potential>::serialize(potential, potential_file);
    JGAP_LOG_INFO("Saved fitted potential to {}", potential_file);

    standardTabulation(potential, output_prefix);

    std::cout << "Execution time: " << formatDuration(elapsedMillisSince(start)) << std::endl;
    return 0;
}
