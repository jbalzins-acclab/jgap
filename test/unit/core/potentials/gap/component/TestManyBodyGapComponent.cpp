#include <cmath>
#include <gtest/gtest.h>
#include "jgap/core/atomic/Atoms.hpp"
#include "jgap/core/atomic/neighbours/NeighbourLists.hpp"
#include "jgap/core/kernels/SquaredExpKernel.hpp"
#include "jgap/core/potentials/gap/component/ManyBodyGapComponent.hpp"

#include "jgap/core/tabulation/TabulationData.hpp"
#include "jgap/core/transform/manybody/TwoBodySum.hpp"
#include "jgap/core/transform/nbody/2b/eam/PolycutoffPairFunction.hpp"

using namespace jgap;

TEST(TestManyBodyGapComponent, RealPolycutoff) {
    Atoms atoms({{0, 0, 0}, {2, 0, 0}}, {Species("Fe"), Species("Ni")});
    auto nl = NeighbourLists(atoms, 5.0);
    auto trans = PolycutoffPairFunction(4.0, 1.0, 2.0);

    auto eam_aggregator = TwoBodySum<1>("Fe");
    eam_aggregator.extend({"Fe", "Ni"}, trans);

    auto kernel = SquaredExpKernel<1, 0>(1.0, std::array<Real, 1>{1.0});
    std::vector<Descriptor<1>> sparse_points = {{1.5}};

    auto component = ManyBodyGapComponent<1, SquaredExpKernel<1, 0>>(eam_aggregator, kernel, sparse_points);

    auto result_opt = component.covariate(nl);
    ASSERT_TRUE(result_opt.has_value());
    auto& quantities = result_opt.value();

    // --- Analytical Calculation ---
    // Central atom: Fe(0). Neighbor: Ni(1). Distance r = 2.0.
    // PolycutoffPairFunction evaluation:
    // chi = (2.0 - 1.0) / (4.0 - 1.0) = 1/3
    // val = 2.0 * (1 - (1/27) * (6(1/9) - 15(1/3) + 10)) = 128/81
    Real expected_val = 128.0 / 81.0;

    // deriv = 2.0 * (1/3) * (1/9) * (-30(1/9) + 60(1/3) - 30) = -80/81
    Real expected_deriv = -80.0 / 81.0;

    // Kernel Value:
    // K = exp(-0.5 * (1.5 - 128/81)^2)
    Real exp_arg = std::pow(1.5 - expected_val, 2.0);
    Real expected_K = std::exp(-0.5 * exp_arg);
    EXPECT_NEAR(quantities.energy(0), expected_K, 1e-9);

    // Kernel Gradient:
    // gradK = K * (1.5 - 128/81)
    Real gradK = expected_K * (1.5 - expected_val);

    // Forces: F_i = -dE/dx_i = (dE/dr) * direction_i
    // dE/dr = (dq/dr) * (dK/dq) = expected_deriv * gradK
    // u = direction from 0 to 1 = (1, 0, 0)
    // F_0 = dE/dr * u. F_1 = -dE/dr * u
    Vector3 expected_force_Fe{expected_deriv * gradK, 0.0, 0.0};
    Vector3 expected_force_Ni{-expected_deriv * gradK, 0.0, 0.0};

    EXPECT_NEAR(quantities.force(0, 0).x, expected_force_Fe.x, 1e-9);
    EXPECT_NEAR(quantities.force(0, 0).y, 0.0, 1e-9);
    EXPECT_NEAR(quantities.force(0, 0).z, 0.0, 1e-9);

    EXPECT_NEAR(quantities.force(0, 1).x, expected_force_Ni.x, 1e-9);
    EXPECT_NEAR(quantities.force(0, 1).y, 0.0, 1e-9);
    EXPECT_NEAR(quantities.force(0, 1).z, 0.0, 1e-9);

    // Virials: V = -F \otimes r
    Real expected_virial_xx = -2.0 * expected_deriv * gradK;
    EXPECT_NEAR(quantities.virials(0).xx, expected_virial_xx, 1e-9);
}

TEST(TestManyBodyGapComponent, TabulationRealPolycutoff) {
    auto trans = PolycutoffPairFunction(4.0, 1.0, 2.0);
    auto eam_aggregator = TwoBodySum<1>("Fe");
    eam_aggregator.extend({"Fe", "Ni"}, trans);

    auto kernel = SquaredExpKernel<1, 0>(1.0, std::array<Real, 1>{1.0});
    std::vector<Descriptor<1>> sparse_points = {{1.5}};
    std::vector<Real> coeffs = {3.0};

    auto component = ManyBodyGapComponent<1, SquaredExpKernel<1, 0>>(eam_aggregator, kernel, sparse_points);

    TabulationParams params;
    params.n_grid_2b = 10;
    params.max_cutoffs = component.getCutoffs();
    params.max_eam_density = 5.0;
    TabulationData tables(params);

    ASSERT_THROW(component.tabulate(tables), std::runtime_error);

    component.setCoefficients(coeffs);
    component.tabulate(tables);

    auto& eam_grids = tables.eam_grids_vec.back();
    auto& value_grid = eam_grids.value_grid;

    for (size_t i = 0; i < value_grid.data_flat.size(); ++i) {
        Real rho = value_grid.getCoord({i})[0];
        Real expected_K = std::exp(-0.5 * std::pow(1.5 - rho, 2));
        EXPECT_NEAR(value_grid.data_flat[i], coeffs[0] * expected_K, 1e-9);
    }
}
