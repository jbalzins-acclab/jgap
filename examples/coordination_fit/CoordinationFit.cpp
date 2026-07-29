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

    return createCoordinationComponents<Dim, SquaredExpKernel<Dim, 0>>(
        coord_trans, kernel_coord, sparsifier_coord, training_data
    );
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

    StandardGapParams params{.seed = 120};
    if (argc == 4) {
        params.zbl_dataset_file = argv[3];
    }
    params.eam_pf = FSGenPairFunction(4.5_r, 3.0_r);
    params.eam_mode = EamMode::Blind;
    params.eam_n_sparse = 20;
    params.n_sparse2 = 20;
    params.n_sparse3 = 500;
    params.regularization_rules = SimpleRegularizationRules();

    GapPotential potential;

    // 2-Body Components
    if (params.n_sparse2 > 0) {
        auto trans2 = PairDistanceTransformation(CosCutoff(params.cutoff2, params.cutoff2_width));
        auto kernel2 = SquaredExpKernel<1, 1>(10.0_r, {1.0_r});
        auto sparsifier2 = HistogramUniformSparsifier<2>(params.seed, params.n_sparse2, std::array{true, false});
        potential.addComponents(
            createTwoBodyComponents<2, SquaredExpKernel<1, 1>>(training_data, trans2, kernel2, sparsifier2)
        );
    }

    // Coordination Components
    // BCC
    potential.addComponents(
        makeCoordinationComponents<4>(
            training_data, params.seed, 500,
            std::array<std::pair<Real, Real>, 4>{{
                {2.4_r, 2.6_r},
                {2.75_r, 2.95_r},
                {3.95_r, 4.15_r},
                {4.65_r, 4.85_r},
            }}
        )
    );

    // FCC
    potential.addComponents(
        makeCoordinationComponents<4>(
            training_data, params.seed, 500,
            std::array<std::pair<Real, Real>, 4>{{
                {2.45_r, 2.75_r},
                {3.5_r, 3.8_r},
                {4.3_r, 4.6_r},
                {5.0_r, 5.3_r},
            }}
        )
    );

    // EAM Components
    if (params.eam_n_sparse > 0) {
        auto kernel_eam = SquaredExpKernel<1, 0>(1.0_r, {1.0_r});
        auto sparsifier_eam = HistogramUniformSparsifier<1>(params.seed, params.eam_n_sparse);
        potential.addComponents(
            createEamComponents<SquaredExpKernel<1, 0>>(
                params.eam_pf, kernel_eam, sparsifier_eam, training_data, params.eam_mode
            )
        );
    }

    // 3-Body Components
    if (params.n_sparse3 > 0) {
        auto trans3 = Angle3bTransformation(CosCutoff(params.cutoff3, params.cutoff3_width));
        auto kernel3 = SquaredExpKernel<3, 1>(1.0_r, {1.0_r, 1.0_r, 1.0_r});
        auto sparsifier3 =
            HistogramUniformSparsifier<4>(params.seed, params.n_sparse3, std::array{true, true, true, false});
        potential.addComponents(
            createThreeBodyComponents<4, SquaredExpKernel<3, 1>>(training_data, trans3, kernel3, sparsifier3)
        );
    }

    IsolatedAtomPotential isolated_atom_pot{training_data};
    ZblPotential zbp_pot =
        params.zbl_dataset_file ? ZblPotential{*params.zbl_dataset_file, training_data} : ZblPotential{training_data};

    CompositePotential external{{
        {"isolated", isolated_atom_pot},
        {"zbl", zbp_pot},
    }};

    potential.optional_external_potential = external;

    // Fit
    QRGapFit fitter;
    fitter.fit(potential, training_data, *params.regularization_rules);

    const std::string potential_file = output_prefix + ".jgap.h5";
    SerializationRegistry<Potential>::serialize(potential, potential_file);
    JGAP_LOG_INFO("Saved fitted potential to {}", potential_file);

    std::cout << "Execution time: " << formatDuration(elapsedMillisSince(start)) << std::endl;
    return 0;
}
