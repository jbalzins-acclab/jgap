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
#include <set>
#include <string>

#include "core/atomic/Atoms.hpp"
#include "core/atomic/neighbours/NeighbourList.hpp"
#include "core/atomic/species/Species.hpp"
#include "core/atomic/species/SpeciesSet.hpp"
#include "core/cutoff/CosCutoff.hpp"
#include "core/kernels/SquaredExpKernel.hpp"
#include "core/sparsification/HistogramUniformSparsifier.hpp"
#include "core/transform/2b/TwoBodyTransformation.hpp"
#include "core/transform/3b/Angle3bTransformation.hpp"
#include "core/transform/aggregated/TransformationAggregatorImpl.hpp"
#include "core/transform/eam/FSGenPairFunction.hpp"
#include "core/transform/eam/PolycutoffPairFunction.hpp"
#include "core/transform/eam/CoscutoffPairFunction.hpp"
#include "core/potentials/CompositePotential.hpp"
#include "core/potentials/gap/GapPotential.hpp"
#include "core/potentials/gap/component/ManyBodyGapComponent.hpp"
#include "core/potentials/gap/component/NBodyGapComponent.hpp"
#include "core/potentials/isolated/IsolatedAtomPotential.hpp"
#include "core/potentials/zbl/ZblPotential.hpp"
#include "core/fit/gap/QRGapFit.hpp"
#include "core/fit/gap/regularization/SimpleRegularizationRules.hpp"
#include "core/ValuePtr.hpp"
#include "core/fit/gap/regularization/PerConfigTypeRegularizationRules.hpp"
#include "serialization/SerializationRegistry.hpp"
#include "utils/Utils.hpp"
#include "io/log/CurrentLogger.hpp"

using namespace jgap;

namespace {
    // ---- baseline parameters, from test-local/fit-params.json ----
    constexpr size_t SEED = 120;

    constexpr Real EAM_CUTOFF = 5.0;          // json eam pair_function cutoff
    constexpr Real EAM_ENERGY_SCALE = 1.0;    // json eam kernel_setup energy_scale
    constexpr size_t EAM_N_SPARSE = 20;       // json eam kernel_setup n_sparse

    constexpr Real CUTOFF_3B = 4.0;           // json 3b cutoff
    constexpr Real WIDTH_3B = 0.6;            // json 3b cutoff_transition_width
    constexpr Real ENERGY_SCALE_3B = 1.0;     // json 3b energy_scale (base, before the Ni adjustment)
    constexpr size_t N_SPARSE_3B = 500;       // json 3b n_sparse (per triplet — reduce for a quick run)

    constexpr Real CUTOFF_2B = 5.0;           // json 2b cutoff
    constexpr Real WIDTH_2B = 0.5;            // json 2b cutoff_transition_width
    constexpr Real ENERGY_SCALE_2B = 10.0;    // json 2b energy_scale
    constexpr size_t N_SPARSE_2B = 20;        // json 2b n_sparse

    bool isNi(const Species& s) { return s.symbol() == "Ni"; }

    /// Picks a different EAM pair (density) function per element pair, to exercise all three kinds.
    ValuePtr<NBodyTransformation<1, 2>> makeEamPairFunction(const Species& central, const Species& contributor) {
        const int n_ni = (isNi(central) ? 1 : 0) + (isNi(contributor) ? 1 : 0);
        if (n_ni == 0) {
            return FSGenPairFunction(EAM_CUTOFF, /*degree=*/3.0);   // Fe-Fe
        }
        if (n_ni == 2) {
            return CoscutoffPairFunction(EAM_CUTOFF, /*r_min=*/0.0); // Ni-Ni
        }
        return PolycutoffPairFunction(EAM_CUTOFF, /*r_min=*/0.0);     // mixed Fe-Ni
    }
}

int main(int argc, char** argv) {
    CurrentLogger::initDefault({.stdout_log_debug = true});

    const std::string training_file = argc > 1 ? argv[1] : "resources/xyz-samples/feni-train.xyz";
    const std::string output_prefix = argc > 2 ? argv[2] : "feni-custom";

    const auto start = std::chrono::steady_clock::now();

    JGAP_LOG_INFO("Custom FeNi fit on {}", training_file);
    const auto training_data = Atoms::readAtoms(training_file);

    // Elements present in the training data.
    std::set<Species> elements;
    for (const auto& atoms : training_data) {
        for (const auto& s : atoms.lookupSpecies()) {
            elements.insert(s);
        }
    }

    GapPotential potential;

    // ===== EAM: one component per central species, with a per-pair density function =====
    const auto eam_kernel = SquaredExpKernel<1, 0>(EAM_ENERGY_SCALE, {1.0});
    const HistogramUniformSparsifier<1> eam_sparsifier(SEED, EAM_N_SPARSE);
    for (const Species& central : elements) {
        auto aggregator = std::make_unique<TransformationAggregatorImpl<1, 2>>(central);
        for (const Species& contributor : elements) {
            aggregator->extend({central, contributor}, makeEamPairFunction(central, contributor));
        }
        potential.addComponent(ManyBodyGapComponent(
            ValuePtr<TransformationAggregator<1>>(std::move(aggregator)),
            eam_kernel, eam_sparsifier, training_data));
    }

    // ===== 3-body: per-triplet energy scale (grows slightly with the Ni count) =====
    const ValuePtr<NBodyTransformation<4, 3>> trans3 = Angle3bTransformation(CosCutoff(CUTOFF_3B, WIDTH_3B));
    const HistogramUniformSparsifier<4> sparsifier3(SEED, N_SPARSE_3B, std::array{true, true, true, false});

    std::set<SpeciesSet<3, HasCentralAtom>> triplets;
    for (const auto& atoms : training_data) {
        NeighbourList nl(atoms, trans3->getCutoffs().maxOverall());
        auto sets = nl.getSpeciesSets<3, HasCentralAtom>();
        triplets.insert(sets.begin(), sets.end());
    }
    for (const auto& triplet : triplets) {
        const int n_ni = (isNi(triplet.getRoot()) ? 1 : 0)
                       + (isNi(triplet.getNodes()[0]) ? 1 : 0)
                       + (isNi(triplet.getNodes()[1]) ? 1 : 0);
        const Real energy_scale = ENERGY_SCALE_3B * (1.0 + 0.1 * static_cast<Real>(n_ni)); // mild Ni boost
        const auto kernel3 = SquaredExpKernel<3, 1>(energy_scale, {1.0, 1.0, 1.0});

        potential.addComponent(NBodyGapComponent(triplet, trans3, kernel3, sparsifier3, training_data));
    }

    // ===== 2-body: one shared setup (as in the standard fit) =====
    const ValuePtr<NBodyTransformation<2, 2>> trans2 = TwoBodyTransformation(CosCutoff(CUTOFF_2B, WIDTH_2B));
    const auto kernel2 = SquaredExpKernel<1, 1>(ENERGY_SCALE_2B, {1.0});
    const HistogramUniformSparsifier<2> sparsifier2(SEED, N_SPARSE_2B, std::array{true, false});
    potential.addComponents(NBodyGapComponent<2, 2, FullSymmetry, SquaredExpKernel<1, 1>>::createComponents(
        training_data, trans2, kernel2, sparsifier2));

    // ===== external: isolated-atom energies + ZBL repulsion =====
    potential.optional_external_potential = CompositePotential{{
        {"isolated", IsolatedAtomPotential{training_data}},
        {"zbl", ZblPotential{training_data}},
    }};

    // ===== setup regularization =====
    const PerConfigTypeRegularizationRules regularization(
        ConfigSigmas(0.002, 0.1, 0.2),
        "isolated_atom:0.0001:0.04:0.04:0.0:liquid:0.01:0.5:2.0:0.0:dimer:0.01:0.5:2.0:0.0:"
                         "short_range:0.01:0.5:2.0:0.0:liquid_surface_100:0.01:0.5:2.0:0.0:"
                         "liquid_surface_110:0.01:0.5:2.0:0.0:liquid_surface_111:0.01:0.5:2.0:0.0:"
                         "gamma_surface:0.002:0.08:0.5:0.0:liquid_high:0.02:0.8:5.0:0.0:"
                         "binary_alloy_melting:0.01:0.5:2.0:0.0:binary_alloy_short_range:0.01:0.5:2.0:0.0"
    );

    // ===== fit =====
    QRGapFit fitter;
    fitter.fit(potential, training_data, regularization);

    const std::string potential_file = output_prefix + ".jgap.h5";
    SerializationRegistry<Potential>::serialize(potential, potential_file);
    JGAP_LOG_INFO("Saved fitted potential to {}", potential_file);

    std::cout << "Execution time: " << formatDuration(elapsedMillisSince(start)) << std::endl;
    return 0;
}
