// Example: fit a standard 2b+3b+EAM GAP on a training set, serialize it, and tabulate it.
//
// Usage: standard_fit <training.xyz> <output_prefix>
//   writes <output_prefix>.jgap.h5 (the serialized potential) and, via TabGapIO,
//   <output_prefix>.tabgap.h5 + <output_prefix>.eam.fs file(s).

#include <chrono>
#include <iostream>
#include <string>

#include "jgap/jgap.hpp"

using namespace jgap;

int main(int argc, char** argv) {
    CurrentLogger::initDefault({.stdout_log_debug = true});

    if (argc < 3 || argc > 4) {
        std::cerr << "Usage: " << argv[0] << " <training.xyz> <output_prefix> [zbl_dataset_file]\n";
        return 1;
    }
    const std::string training_file = argv[1];
    const std::string output_prefix = argv[2];

    const auto start = std::chrono::steady_clock::now();

    JGAP_LOG_INFO("Fitting on {}", training_file);
    auto training_data = Atoms::readAtoms(training_file);

    StandardGapParams params{.seed = 120};
    if (argc == 4) {
        params.zbl_dataset_file = argv[3]; // otherwise the built-in ZBL dataset is used
    }
    params.eam_pf = FSGenPairFunction(4.5, 3.0);
    params.eam_mode = EamMode::Blind;
    params.n_sparse3 = 500;
    params.regularization_rules = SimpleRegularizationRules();

    auto potential = standardGapFit(training_data, params);

    // Serialize the fitted potential; the registry stamps the node with the concrete type so it can be
    // read back (see read_and_predict) without knowing what kind of potential it is.
    const std::string potential_file = output_prefix + ".jgap.h5";
    SerializationRegistry<Potential>::serialize(potential, potential_file);
    JGAP_LOG_INFO("Saved fitted potential to {}", potential_file);

    TabulationData tabulation_data = potential.tabulate(
        {.max_cutoffs = potential.getCutoffs(), .max_eam_density = 10.0, .n_grid_2b = 5000, .n_grid_3b = {80, 80, 80}}
    );
    TabGapPotential tabgap{tabulation_data};

    const Filenames tabgap_files = TabGapIO::write(tabgap, output_prefix);
    JGAP_LOG_INFO("Saved tabGAP to: {}", vectorToString(tabgap_files));

    std::cout << "Execution time: " << formatDuration(elapsedMillisSince(start)) << std::endl;
    return 0;
}
