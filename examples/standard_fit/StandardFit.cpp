// Example: fit a standard 2b+3b+EAM GAP on a training set, serialize it, and tabulate it.
//
// Usage: standard_fit <training.xyz> <output_prefix> [screened_coulomb_dataset_file]
//   writes <output_prefix>.jgap.h5 (the serialized potential) and, via standardTabulation,
//   <output_prefix>.tabgap.h5 + <output_prefix>.eam.fs file(s).

#include <chrono>
#include <iostream>
#include <string>

#include "jgap/core/fit/gap/regularization/PerConfigTypeRegularizationRules.hpp"
#include "jgap/jgap.hpp"

using namespace jgap;
using namespace jgap::utils;

int main(int argc, char** argv) {
    CurrentLogger::initDefault({.stdout_log_debug = true});

    if (argc < 3 || argc > 4) {
        std::cerr << "Usage: " << argv[0] << " <training.xyz> <output_prefix> [screened_coulomb_dataset_file]\n";
        return 1;
    }
    const std::string training_file = argv[1];
    const std::string output_prefix = argv[2];

    const auto total_start = std::chrono::steady_clock::now();

    JGAP_LOG_INFO("Fitting on {}", training_file);
    auto training_data = Atoms::readAtoms(training_file);

    StandardGapParams params{.seed = 120};
    if (argc == 4) {
        params.screened_coulomb_dataset_file = argv[3]; // otherwise the built-in screening dataset is used
    }
    params.eam_pair_function = EamPairFunctionType::FSGen3;
    params.eam_mode = EamMode::Blind;
    params.n_sparse3 = 500;

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
