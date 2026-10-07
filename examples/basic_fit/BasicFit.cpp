// Example: fit a 2b+3b+EAM GAP with SquaredExpKernel assembled manually without utils,
// serialize it, and tabulate it with a fast 12x12x12 3-body grid.
//
// Usage: basic_fit <training.xyz> <output_prefix> [screened_coulomb_dataset_file] [--ram-limit <gb>]
//   writes <output_prefix>.jgap.h5 (the serialized potential) and, via standardTabulation,
//   <output_prefix>.tabgap.h5 + <output_prefix>.eam.fs file(s).

#include <chrono>
#include <iostream>
#include <optional>
#include <set>
#include <string>
#include <vector>

#include "jgap/core/atomic/species/composition/Species2Sorted.hpp"
#include "jgap/core/atomic/species/composition/Species3AtomicSorted.hpp"
#include "jgap/core/fit/gap/regularization/PerConfigTypeRegularizationRules.hpp"
#include "jgap/impl/transform/manybody/TwoBodySum.hpp"
#include "jgap/impl/fit/gap/ElementIncrementalQRGapFit.hpp"
#include "jgap/jgap.hpp"

using namespace jgap;
using namespace jgap::utils;

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
            std::cout << "Usage: " << argv[0] << " <training.xyz> <output_prefix> [screened_coulomb_dataset_file] [--ram-limit <gb>]\n"
                      << "  --ram-limit <gb>   RAM limit in GB for ElementIncrementalQRGapFit (default: 2.0)\n";
            return 0;
        } else {
            positional_args.push_back(arg);
        }
    }

    if (positional_args.size() < 2 || positional_args.size() > 3) {
        std::cerr << "Usage: " << argv[0] << " <training.xyz> <output_prefix> [screened_coulomb_dataset_file] [--ram-limit <gb>]\n";
        return 1;
    }
    const std::string training_file = positional_args[0];
    const std::string output_prefix = positional_args[1];
    const std::optional<std::string> screened_coulomb_dataset_file =
        (positional_args.size() == 3) ? std::optional<std::string>{positional_args[2]} : std::nullopt;

    const auto start = std::chrono::steady_clock::now();

    JGAP_LOG_INFO("Manual assembly fit on {} using ElementIncrementalQRGapFit (RAM limit: {} GB)", training_file, ram_limit_gb);
    auto training_data = Atoms::readAtoms(training_file);

    const size_t seed = 120;

    // 2-body parameters
    const double cutoff2 = 4.5;
    const double cutoff2_width = 1.0;
    const size_t n_sparse2 = 20;

    // EAM parameters
    const auto eam_pf = FSGenPairFunction(cutoff2, 3.0);
    const EamMode eam_mode = EamMode::Blind;
    const size_t eam_n_sparse = 20;

    // 3-body parameters
    const double cutoff3 = 3.7;
    const double cutoff3_width = 0.6;
    const size_t n_sparse3 = 500;

    // Discover species from training data
    std::set<Species> elements;
    for (const auto& atoms : training_data) {
        for (const auto& s : atoms.getSpecies()) {
            elements.insert(s);
        }
    }

    GapPotential potential;

    // ====================================================================================
    // 1. 2-Body Components with SquaredExpKernel (assembled directly without utils)
    // ====================================================================================
    if (n_sparse2 > 0) {
        auto trans2 = PairDistanceTransformation(CosCutoff(cutoff2, cutoff2_width));
        auto kernel2 = SquaredExpKernel<1, 1>(10.0, {1.0});
        auto sparsifier2 = HistogramUniformSparsifier<2>(seed, n_sparse2, std::array{true, false});

        std::set<Species2Sorted> pairs;
        for (const auto& atoms : training_data) {
            NeighbourLists nl(atoms, cutoff2);
            auto sets = Species2Sorted::getAll(nl);
            pairs.insert(sets.begin(), sets.end());
        }
        for (const auto& pair : pairs) {
            potential.addComponent(TwoBodyGapComponent<2, SquaredExpKernel<1, 1>>(
                pair, trans2, kernel2, sparsifier2, training_data
            ));
        }
    }

    // ====================================================================================
    // 2. ManyBodyGapComponent with EAM Pair Function and SquaredExpKernel (without utils)
    // ====================================================================================
    if (eam_n_sparse > 0) {
        auto kernel_eam = SquaredExpKernel<1, 0>(1.0, {1.0});
        auto sparsifier_eam = HistogramUniformSparsifier<1>(
            seed, eam_n_sparse, std::nullopt, std::nullopt, Descriptor<1>{0.05}
        );

        for (const auto& central : elements) {
            auto aggregator = std::make_unique<TwoBodySum<1>>(central);
            for (const auto& contributor : elements) {
                aggregator->extend({central, contributor}, eam_pf);
            }
            potential.addComponent(ManyBodyGapComponent<1, SquaredExpKernel<1, 0>>(
                ValuePtr<NBodyAggregator<1>>(std::move(aggregator)),
                kernel_eam,
                sparsifier_eam,
                training_data
            ));
        }
    }

    // ====================================================================================
    // 3. 3-Body Components with SquaredExpKernel (assembled directly without utils)
    // ====================================================================================
    if (n_sparse3 > 0) {
        auto trans3 = Angle3bTransformation(CosCutoff(cutoff3, cutoff3_width));
        auto kernel3 = SquaredExpKernel<3, 1>(1.0, {1.0, 1.0, 1.0});
        auto sparsifier3 = HistogramUniformSparsifier<4>(seed, n_sparse3, std::array{true, true, true, false});

        std::set<Species3AtomicSorted> triplets;
        for (const auto& atoms : training_data) {
            NeighbourLists nl(atoms, cutoff3);
            auto sets = Species3AtomicSorted::getAll(nl);
            triplets.insert(sets.begin(), sets.end());
        }
        for (const auto& triplet : triplets) {
            potential.addComponent(ThreeBodyGapComponent<4, SquaredExpKernel<3, 1>>(
                triplet, trans3, kernel3, sparsifier3, training_data
            ));
        }
    }

    // External potentials
    IsolatedAtomPotential isolated_atom_pot{training_data};
    ScreenedCoulombPotential sc_pot = screened_coulomb_dataset_file
                                          ? ScreenedCoulombPotential{*screened_coulomb_dataset_file, training_data}
                                          : ScreenedCoulombPotential{training_data};

    CompositePotential external{{
        {"isolated", isolated_atom_pot},
        {"screened_coulomb", sc_pot},
    }};

    potential.optional_external_potential = external;

    // Fit with ElementIncrementalQRGapFit
    ElementIncrementalQRGapFit fitter(ram_limit_gb);
    PerConfigTypeRegularizationRules regularization_rules{PerConfigTypeSigmas(0.001, 0.05, 0.1, 0.02)};
    auto sigmas = regularization_rules.determineForAll(training_data);
    fitter.fit(potential, training_data, sigmas);

    const std::string potential_file = output_prefix + ".jgap.h5";
    SerializationRegistry<Potential>::serialize(potential, potential_file);
    JGAP_LOG_INFO("Saved fitted potential to {}", potential_file);

    // Tabulation with custom 12x12x12 3B grid
    StandardTabulationParams tab_params;
    tab_params.n_grid_3b = {12, 12, 12};
    standardTabulation(potential, output_prefix, tab_params);

    std::cout << "Execution time: " << formatDuration(elapsedMillisSince(start)) << std::endl;
    return 0;
}
