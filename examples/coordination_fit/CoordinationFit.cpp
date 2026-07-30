#include <chrono>
#include <iostream>
#include <string>

#include "jgap/experimental/transform/nbody/2b/CoordinationTransformation.hpp"
#include "jgap/jgap.hpp"
#include "jgap/utils/gap/GapComponentUtils.hpp"

using namespace jgap;

template<size_t Dim>
auto makeCoordinationComponents(
    const std::vector<Atoms>& training_data, uint64_t seed, size_t n_sparse,
    const std::array<std::pair<Real, Real>, Dim>& ranges
) {
    auto coord_trans = ValuePtr<CoordinationTransformation<Dim>>(CoordinationTransformation<Dim>(ranges));

    std::array<Real, Dim> length_scales;
    length_scales.fill(1.0_r);
    auto kernel_coord = SquaredExpKernel<Dim, 0>(1.0_r, length_scales);

    std::array<bool, Dim> use_histogram;
    use_histogram.fill(true);
    auto sparsifier_coord = HistogramUniformSparsifier<Dim>(seed, n_sparse, use_histogram);

    return createCoordinationComponents(coord_trans, kernel_coord, sparsifier_coord, training_data);
}

auto makeMeamComponents(const std::vector<Atoms>& training_data, uint64_t seed, size_t n_sparse, Real cutoff) {
    auto meam_trans = ValuePtr<MeamTransformation>(MeamTransformation(CosCutoff(cutoff, 0.5_r)));

    auto kernel_coord = SquaredExpKernel<3, 0>(1.0_r, {1.0_r, 1.0_r, 1.0_r});

    auto sparsifier_coord = HistogramUniformSparsifier<3>(seed, n_sparse, std::array{true, true, true});

    return createMeamComponents(meam_trans, kernel_coord, sparsifier_coord, training_data);
}

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

    std::optional<std::string> zbl_dataset_file;
    if (argc == 4) {
        zbl_dataset_file = argv[3];
    }

    GapPotential potential;

    // 2-Body Components
    if (true) {
        auto trans2 = PairDistanceTransformation(CosCutoff(4.5_r, 1.0_r));
        auto kernel2 = SquaredExpKernel<1, 1>(10.0_r, {1.0_r});
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
                training_data, 120, 500,
                std::array<std::pair<Real, Real>, 3>{{
                    {2.4_r, 2.6_r},
                    {2.75_r, 2.95_r},
                    {3.95_r, 4.15_r},
                    //{4.65_r, 4.85_r},
                }}
            )
        );
    }
    if (false) {
        // FCC
        potential.addComponents(
            makeCoordinationComponents<3>(
                training_data, 120, 500,
                std::array<std::pair<Real, Real>, 3>{{
                    {2.45_r, 2.75_r},
                    {3.5_r, 3.8_r},
                    {4.3_r, 4.6_r},
                    //{5.0_r, 5.3_r},
                }}
            )
        );
    }

    // MEAM Components
    if (true) {
        potential.addComponents(makeMeamComponents(training_data, 120, 500, 3.7_r));
    }

    // EAM Components
    if (true) {
        auto kernel_eam = SquaredExpKernel<1, 0>(1.0_r, {1.0_r});
        auto sparsifier_eam = HistogramUniformSparsifier<1>(120, 20);
        potential.addComponents(
            createEamComponents<SquaredExpKernel<1, 0>>(
                FSGenPairFunction(4.5_r, 3.0_r), kernel_eam, sparsifier_eam, training_data, EamMode::Blind
            )
        );
    }

    // 3-Body Components
    if (true) {
        auto trans3 = Angle3bTransformation(CosCutoff(3.7_r, 0.6_r));
        auto kernel3 = SquaredExpKernel<3, 1>(1.0_r, {1.0_r, 1.0_r, 1.0_r});
        auto sparsifier3 = HistogramUniformSparsifier<4>(120, 500, std::array{true, true, true, false});
        potential.addComponents(
            createThreeBodyComponents<4, SquaredExpKernel<3, 1>>(training_data, trans3, kernel3, sparsifier3)
        );
    }

    IsolatedAtomPotential isolated_atom_pot{training_data};
    ZblPotential zbp_pot =
        zbl_dataset_file ? ZblPotential{*zbl_dataset_file, training_data} : ZblPotential{training_data};

    CompositePotential external{{
        {"isolated", isolated_atom_pot},
        {"zbl", zbp_pot},
    }};

    potential.optional_external_potential = external;

    // Fit
    QRGapFit fitter;
    fitter.fit(potential, training_data, SimpleRegularizationRules());

    const std::string potential_file = output_prefix + ".jgap.h5";
    SerializationRegistry<Potential>::serialize(potential, potential_file);
    JGAP_LOG_INFO("Saved fitted potential to {}", potential_file);

    std::cout << "Execution time: " << formatDuration(elapsedMillisSince(start)) << std::endl;
    return 0;
}
