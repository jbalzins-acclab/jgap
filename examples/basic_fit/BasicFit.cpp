// Example: fit a 2b+3b+EAM GAP using QuadraticBumpKernel on a training set, serialize it, and tabulate it.
//
// Usage: basic_fit <training.xyz> <output_prefix> [zbl_dataset_file]
//   writes <output_prefix>.jgap.h5 (the serialized potential) and, via TabGapIO,
//   <output_prefix>.tabgap.h5 + <output_prefix>.eam.fs file(s).

#include <chrono>
#include <iostream>
#include <optional>
#include <string>

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
    const std::optional<std::string> screened_coulomb_dataset_file =
        (argc == 4) ? std::optional<std::string>{argv[3]} : std::nullopt;

    const auto start = std::chrono::steady_clock::now();

    JGAP_LOG_INFO("Fitting on {}", training_file);
    auto training_data = Atoms::readAtoms(training_file);

    const size_t seed = 120;

    // 2-body parameters
    const Real cutoff2 = 4.5;
    const Real cutoff2_width = 1.0;
    const size_t n_sparse2 = 20;

    // EAM parameters
    const auto eam_pf = FSGenPairFunction(4.5, 3.0);
    const EamMode eam_mode = EamMode::Blind;
    const size_t eam_n_sparse = 20;

    // 3-body parameters
    const Real cutoff3 = 3.7;
    const Real cutoff3_width = 0.6;
    const size_t n_sparse3 = 500;

    const SimpleRegularizationRules regularization_rules{};

    GapPotential potential;

    // ====================================================================================
    // 2-Body Components with WendlandKernel
    // ====================================================================================
    if (n_sparse2 > 0) {
        auto trans2 = PairDistanceTransformation(CosCutoff(cutoff2, cutoff2_width));
        auto kernel2 = WendlandKernel<1, 1>(10.0_r, {1.0_r});
        auto sparsifier2 = HistogramUniformSparsifier<2>(seed, n_sparse2, std::array{true, false});
        potential.addComponents(
            TwoBodyGapComponent<2, WendlandKernel<1, 1>>::createComponents(training_data, trans2, kernel2, sparsifier2)
        );
    }

    // ====================================================================================
    // ManyBodyGapComponent with EAM Pair Function and WendlandKernel
    // ====================================================================================
    if (eam_n_sparse > 0) {
        auto kernel_eam = WendlandKernel<1, 0>(1.0_r, {1.0_r});
        auto sparsifier_eam = HistogramUniformSparsifier<1>(seed, eam_n_sparse);
        potential.addComponents(
            EamPairFunction::createComponents(eam_pf, kernel_eam, sparsifier_eam, training_data, eam_mode)
        );
    }

    // ====================================================================================
    // 3-Body Components with WendlandKernel
    // ====================================================================================
    if (n_sparse3 > 0) {
        auto trans3 = Angle3bTransformation(CosCutoff(cutoff3, cutoff3_width));
        auto kernel3 = WendlandKernel<3, 1>(1.0_r, {1.0_r, 1.0_r, 1.0_r});
        auto sparsifier3 = HistogramUniformSparsifier<4>(seed, n_sparse3, std::array{true, true, true, false});
        potential.addComponents(
            ThreeBodyGapComponent<4, WendlandKernel<3, 1>>::createComponents(
                training_data, trans3, kernel3, sparsifier3
            )
        );
    }

    IsolatedAtomPotential isolated_atom_pot{training_data};
    ScreenedCoulombPotential sc_pot = screened_coulomb_dataset_file
                                          ? ScreenedCoulombPotential{*screened_coulomb_dataset_file, training_data}
                                          : ScreenedCoulombPotential{training_data};

    CompositePotential external{{
        {"isolated", isolated_atom_pot},
        {"screened_coulomb", sc_pot},
    }};

    potential.optional_external_potential = external;

    // Fit
    IterativeGapFit fitter;
    fitter.fit(potential, training_data, regularization_rules);

    // Serialize the fitted potential; the registry stamps the node with the concrete type so it can be
    // read back without knowing what kind of potential it is.
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
