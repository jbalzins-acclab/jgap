#include <chrono>
#include <iostream>
#include <string>

#include "jgap/experimental/transform/nbody/2b/CoordinationTransformation.hpp"
#include "jgap/jgap.hpp"
#include "jgap/utils/gap/GapComponentUtils.hpp"

using namespace jgap;
using namespace jgap::utils;

template<size_t Dim>
auto makeCoordinationComponents(
    const std::vector<Atoms>& training_data,
    uint64_t seed,
    size_t n_sparse,
    const std::array<std::pair<double, double>, Dim>& ranges
) {
    auto coord_trans = ValuePtr<CoordinationTransformation<Dim>>(CoordinationTransformation<Dim>(ranges));

    std::array<double, Dim> length_scales;
    length_scales.fill(1.0);
    auto kernel_coord = SquaredExpKernel<Dim, 0>(1.0, length_scales);

    std::array<bool, Dim> use_histogram;
    use_histogram.fill(true);
    auto sparsifier_coord = HistogramUniformSparsifier<Dim>(seed, n_sparse, use_histogram);

    return createCoordinationComponents(coord_trans, kernel_coord, sparsifier_coord, training_data);
}

auto makeMeamComponents(const std::vector<Atoms>& training_data, uint64_t seed, size_t n_sparse, double cutoff) {
    auto meam_trans = ValuePtr<MeamTransformation>(MeamTransformation(CosCutoff(cutoff, 0.5)));

    auto kernel_coord = SquaredExpKernel<3, 0>(1.0, {1.0, 1.0, 1.0});

    auto sparsifier_coord = HistogramUniformSparsifier<3>(seed, n_sparse, std::array{true, true, true});

    return createMeamComponents(meam_trans, kernel_coord, sparsifier_coord, training_data);
}

int main(int argc, char** argv) {
    CurrentLogger::initDefault({.stdout_log_debug = true});

    if (argc < 3 || argc > 4) {
        std::cerr << "Usage: " << argv[0] << " <training.xyz> <output_prefix> [screened_coulomb_dataset_file]\n";
        return 1;
    }
    const std::string training_file = argv[1];
    const std::string output_prefix = argv[2];

    const auto start = std::chrono::steady_clock::now();

    JGAP_LOG_INFO("Fitting on {}", training_file);
    auto training_data = Atoms::readAtoms(training_file);

    std::optional<std::string> screened_coulomb_dataset_file;
    if (argc == 4) {
        screened_coulomb_dataset_file = argv[3];
    }

    GapPotential potential;

    // 2-Body Components
    if (true) {
        auto trans2 = PairDistanceTransformation(CosCutoff(4.5, 1.0));
        auto kernel2 = SquaredExpKernel<1, 1>(10.0, {1.0});
        auto sparsifier2 = HistogramUniformSparsifier<2>(120, 20, std::array{true, false});
        potential.addComponents(
            createTwoBodyComponents<2, SquaredExpKernel<1, 1>>(training_data, trans2, kernel2, sparsifier2)
        );
    }

    // Coordination Components
    if (false) {
        // BCC
        potential.addComponents(
            makeCoordinationComponents<3>(
                training_data,
                120,
                500,
                std::array<std::pair<double, double>, 3>{{
                    {2.4, 2.6},
                    {2.75, 2.95},
                    {3.95, 4.15},
                    //{4.65, 4.85},
                }}
            )
        );
    }
    if (false) {
        // FCC
        potential.addComponents(
            makeCoordinationComponents<3>(
                training_data,
                120,
                500,
                std::array<std::pair<double, double>, 3>{{
                    {2.45, 2.75},
                    {3.5, 3.8},
                    {4.3, 4.6},
                    //{5.0, 5.3},
                }}
            )
        );
    }

    // MEAM Components
    if (true) {
        potential.addComponents(makeMeamComponents(training_data, 120, 500, 3.7));
    }

    // EAM Components
    if (false) {
        auto kernel_eam = SquaredExpKernel<1, 0>(1.0, {1.0});
        auto sparsifier_eam = HistogramUniformSparsifier<1>(120, 20);
        potential.addComponents(
            createEamComponents<SquaredExpKernel<1, 0>>(
                FSGenPairFunction(4.5, 3.0), kernel_eam, sparsifier_eam, training_data, EamMode::Blind
            )
        );
    }

    // 3-Body Components
    if (true) {
        auto trans3 = Angle3bTransformation(CosCutoff(3.7, 0.6));
        auto kernel3 = SquaredExpKernel<3, 1>(1.0, {1.0, 1.0, 1.0});
        auto sparsifier3 = HistogramUniformSparsifier<4>(120, 500, std::array{true, true, true, false});
        potential.addComponents(
            createThreeBodyComponents<4, SquaredExpKernel<3, 1>>(training_data, trans3, kernel3, sparsifier3)
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
    QRGapFit fitter;
    auto sigmas = SimpleRegularizationRules().determineForAll(training_data);
    fitter.fit(potential, training_data, sigmas);

    const std::string potential_file = output_prefix + ".jgap.h5";
    SerializationRegistry<Potential>::serialize(potential, potential_file);
    JGAP_LOG_INFO("Saved fitted potential to {}", potential_file);

    std::cout << "Execution time: " << formatDuration(elapsedMillisSince(start)) << std::endl;
    return 0;
}
