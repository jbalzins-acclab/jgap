// Example: fit a standard 2b+3b+EAM GAP on a training set, serialize it, and tabulate it.
//
// Usage: standard_fit <training.xyz> <output_prefix> [screened_coulomb_dataset_file] [--ram-limit <gb>]
//   writes <output_prefix>.jgap.h5 (the serialized potential) and, via standardTabulation,
//   <output_prefix>.tabgap.h5 + <output_prefix>.eam.fs file(s).

#include <chrono>
#include <iostream>
#include <string>
#include <vector>

#include "jgap/core/fit/gap/regularization/PerConfigTypeRegularizationRules.hpp"
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

    const auto total_start = std::chrono::steady_clock::now();

    JGAP_LOG_INFO("Fitting on {} using ElementIncrementalQRGapFit (RAM limit: {} GB)", training_file, ram_limit_gb);
    auto training_data = Atoms::readAtoms(training_file, {.virials = "virial"});

    StandardGapParams params{.seed = 120};
    if (positional_args.size() == 3) {
        params.screened_coulomb_dataset_file = positional_args[2]; // otherwise the built-in screening dataset is used
    }
    if (params.default_eam) {
        params.default_eam->eam_pair_function = EamPairFunctionType::Polycutoff;
        params.default_eam->eam_mode = EamMode::FSsym;
    }
    if (params.default_3b) {
        params.default_3b->n_sparse = 625;
        // Option: use Distances3bTransformation instead of Angle (default):
        // params.default_3b->transformation_type = ThreeBodyTransformationType::Distances;
    }
    // Default non-species parameters can be disabled by erasing them (setting to std::nullopt):
    // params.default_2b = std::nullopt;
    //
    // Species-specific parameters can be added to species_2b, species_eam, or species_3b:
    // params.species_3b.push_back(StandardGap3bParams{
    //     .species = Species3AtomicSorted("Fe|Fe,Ni"),
    //     .transformation_type = ThreeBodyTransformationType::Distances,
    //     .cutoff = 3.7,
    //     .cutoff_width = 0.6,
    //     .n_sparse = 500,
    // });
    params.approx_ram_limit_gb = ram_limit_gb;

    PerConfigTypeRegularizationRules rules{PerConfigTypeSigmas(0.001, 0.05, 0.1, 0.02)};
    auto sigmas = rules.determineForAll(training_data);

    const std::string potential_file = output_prefix + ".jgap.h5";
    const auto fit_start = std::chrono::steady_clock::now();
    standardGapFit(potential_file, training_data, sigmas, params);
    const auto fit_duration = elapsedMillisSince(fit_start);

    // Tabulate fitted potential
    const auto tab_start = std::chrono::steady_clock::now();
    standardTabulation(potential_file, output_prefix);
    const auto tab_duration = elapsedMillisSince(tab_start);

    const auto total_duration = elapsedMillisSince(total_start);

    std::cout << "Fitting execution time:    " << formatDuration(fit_duration) << "\n";
    std::cout << "Tabulation execution time: " << formatDuration(tab_duration) << "\n";
    std::cout << "Total execution time:      " << formatDuration(total_duration) << std::endl;
    return 0;
}
