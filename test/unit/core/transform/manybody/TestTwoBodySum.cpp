#include <gtest/gtest.h>
#include "jgap/core/atomic/Atoms.hpp"
#include "jgap/core/atomic/neighbours/NeighbourLists.hpp"
#include "jgap/core/tabulation/TabulationData.hpp"
#include "jgap/core/transform/manybody/NBodyAggregator.hpp"
#include "jgap/core/transform/manybody/TwoBodySum.hpp"
#include "jgap/core/transform/nbody/2b/eam/EamPairFunction.hpp"
#include "jgap/core/transform/nbody/2b/eam/PolycutoffPairFunction.hpp"

using namespace jgap;

namespace {
    // Mock EAM function that returns a constant value and a distance-dependent derivative
    class MockEamPairFunction final : public EamPairFunction {
    public:
        MockEamPairFunction(Real value_to_return, Real deriv_factor, Real cutoff = 10.0, bool is_rot_inv = true) :
            EamPairFunction(cutoff),
            value_to_return_(value_to_return),
            deriv_factor_(deriv_factor),
            is_rot_inv(is_rot_inv) {}

        Descriptor<1> evaluate(const Cluster2& pair) const override { return {{value_to_return_}}; }

        TwoBodyDescriptor<1> evaluateAndDifferentiate(const Cluster2& pair) const override {
            Real deriv = pair.separation01.magnitude * deriv_factor_;
            return {{value_to_return_}, {deriv * pair.separation01.direction}};
        }

        bool isRotationallyInvariant() const override { return is_rot_inv; }

        MockEamPairFunction* clone() const override { return new MockEamPairFunction(*this); }

    private:
        Real value_to_return_;
        Real deriv_factor_;
        bool is_rot_inv = true;
    };
}

TEST(TestTwoBodySum, SingleTransformation) {
    Atoms atoms(
        {{0, 0, 0}, {1, 1, 0}, {1, 0, 1}, {0, 1, 1}}, {Species("Fe"), Species("Fe"), Species("Fe"), Species("Fe")}
    );
    auto nl = NeighbourLists(atoms, 5.0);

    auto aggregator = TwoBodySum<1>("Fe");
    aggregator.extend({"Fe", "Fe"}, MockEamPairFunction(1.0, 0.1));

    auto aggregated_descriptors = aggregator.aggregate(nl);

    ASSERT_EQ(aggregated_descriptors.values.size(), 4);

    // --- Analytical Calculation for atom 0 ---
    // Neighbors are at (1,1,0), (1,0,1), (0,1,1). All have dist=sqrt(2).
    // Mock deriv = dist * 0.1 = sqrt(2) * 0.1
    Real expected_deriv = std::sqrt(2.0) * 0.1;

    // Force on atom i is F_i = dE/dr * direction. Aggregator stores -dE/dr.
    // The aggregator stores the force contribution for each descriptor dimension.
    // Here, Dim=1, so we only care about forces[...][0].
    // For neighbor (1,1,0): direction = (1/sqrt(2), 1/sqrt(2), 0).
    // Force on atom 0: +direction * deriv = (0.1, 0.1, 0.0)
    // For neighbor (1,0,1): direction = (1/sqrt(2), 0, 1/sqrt(2)).
    // Force on atom 0: +direction * deriv = (0.1, 0.0, 0.1)
    // For neighbor (0,1,1): direction = (0, 1/sqrt(2), 1/sqrt(2)).
    // Force on atom 0: +direction * deriv = (0.0, 0.1, 0.1)
    // Total Force on atom 0: (0.2, 0.2, 0.2)
    Vector3 expected_force_on_0{0.2, 0.2, 0.2};

    // Virials: V += sep.virials * deriv
    // For neighbor (1,1,0): sep.virials.xx=1, xy=1, yy=1. Scaled by deriv.
    // For neighbor (1,0,1): sep.virials.xx=1, xz=1, zz=1. Scaled by deriv.
    // For neighbor (0,1,1): sep.virials.yy=1, yz=1, zz=1. Scaled by deriv.
    // Total Virials: xx=0.2, xy=0.1, xz=0.1, yy=0.2, yz=0.1, zz=0.2 (all scaled by sqrt(2)*0.1)
    // Wait, sep.virials is r_x*u_x, etc.
    // For (1,1,0): r=(1,1,0), u=(1/sqrt(2), 1/sqrt(2), 0).
    // xx=1*1/sqrt(2), xy=1*1/sqrt(2), yy=1*1/sqrt(2).
    // Scaled by deriv: xx=0.1, xy=0.1, yy=0.1
    // Total Virials: xx=0.2, xy=0.1, xz=0.1, yy=0.2, yz=0.1, zz=0.2
    // (^ Wrong sign convention) => x(-1)
    Virials expected_virials_0{-0.2, -0.1, -0.1, -0.2, -0.1, -0.2};

    EXPECT_NEAR(aggregated_descriptors.values[0][0], 3.0, 1e-9);

    const auto& force_0 = aggregated_descriptors.force(0, 0)[0];
    EXPECT_NEAR(force_0.x, expected_force_on_0.x, 1e-9);
    EXPECT_NEAR(force_0.y, expected_force_on_0.y, 1e-9);
    EXPECT_NEAR(force_0.z, expected_force_on_0.z, 1e-9);

    const auto& virial_0 = aggregated_descriptors.virials[0][0];
    EXPECT_NEAR(virial_0.xx, expected_virials_0.xx, 1e-9);
}

TEST(TestTwoBodySum, MultipleTransformations) {
    Atoms atoms(
        {{0, 0, 0}, {3, 0, 0}, {0, 4, 0}, {3, 4, 0}}, {Species("Fe"), Species("Fe"), Species("Ni"), Species("Ni")}
    );
    auto nl = NeighbourLists(atoms, 10.0);

    auto aggregator = TwoBodySum<1>("Fe");
    aggregator.extend({"Fe", "Fe"}, MockEamPairFunction(1.0, 0.1));
    aggregator.extend({"Fe", "Ni"}, MockEamPairFunction(10.0, 0.2));

    auto aggregated_descriptors = aggregator.aggregate(nl);

    ASSERT_EQ(aggregated_descriptors.values.size(), 2);

    // --- Check Fe atom at index 0 (0,0,0) ---
    // Force on atom i is F_i = dE/dr * direction
    // Neighbors:
    // 1. Fe at (3,0,0): r=3.0, u=(1,0,0). deriv = 3.0 * 0.1 = 0.3
    //    Force contrib: +u * deriv = (+0.3, 0.0, 0.0)
    //    Virials: xx = -3.0 * 1.0 * 0.3 = -0.9
    // 2. Ni at (0,4,0): r=4.0, u=(0,1,0). deriv = 4.0 * 0.2 = 0.8
    //    Force contrib: +u * deriv = (0.0, +0.8, 0.0)
    //    Virials: yy = -4.0 * 1.0 * 0.8 = -3.2
    // 3. Ni at (3,4,0): r=5.0, u=(0.6, 0.8, 0). deriv = 5.0 * 0.2 = 1.0
    //    Force contrib: +u * deriv = (+0.6, +0.8, 0.0)
    //    Virials: xx = -3.0 * 0.6 * 1.0 = -1.8
    //             xy = -3.0 * 0.8 * 1.0 = -2.4
    //             yy = -4.0 * 0.8 * 1.0 = -3.2
    //
    // Total Force: (+0.9, +1.6, 0.0)
    // Total Virials: xx = -2.7, xy = -2.4, xz = 0.0, yy = -6.4, yz = 0.0, zz = 0.0

    EXPECT_NEAR(aggregated_descriptors.values[0][0], 21.0, 1e-9);

    const auto& calc_force0 = aggregated_descriptors.force(0, 0)[0];
    EXPECT_NEAR(calc_force0.x, 0.9, 1e-9);
    EXPECT_NEAR(calc_force0.y, 1.6, 1e-9);
    EXPECT_NEAR(calc_force0.z, 0.0, 1e-9);

    const auto& calc_virials0 = aggregated_descriptors.virials[0][0];
    EXPECT_NEAR(calc_virials0.xx, -2.7, 1e-9);
    EXPECT_NEAR(calc_virials0.xy, -2.4, 1e-9);
    EXPECT_NEAR(calc_virials0.xz, 0.0, 1e-9);
    EXPECT_NEAR(calc_virials0.yy, -6.4, 1e-9);
    EXPECT_NEAR(calc_virials0.yz, 0.0, 1e-9);
    EXPECT_NEAR(calc_virials0.zz, 0.0, 1e-9);
}

TEST(TestTwoBodySum, GetCutoffs) {
    auto aggregator = TwoBodySum<1>("Fe");
    aggregator.extend({"Fe", "Fe"}, MockEamPairFunction(1.0, 0.1, 5.0));
    aggregator.extend({"Fe", "Ni"}, MockEamPairFunction(10.0, 0.2, 8.0));

    Cutoffs cutoffs = aggregator.getCutoffs();
    EXPECT_NEAR(cutoffs.forDim(2), 8.0, 1e-9);
}

TEST(TestTwoBodySum, Tabulation) {
    auto aggregator = TwoBodySum<1>("Fe");
    aggregator.extend({"Fe", "Fe"}, MockEamPairFunction(1.0, 0.1, 5.0));
    aggregator.extend({"Fe", "Ni"}, MockEamPairFunction(10.0, 0.2, 8.0));

    TabulationParams params;
    params.n_grid_2b = 10;
    params.max_cutoffs = aggregator.getCutoffs();
    TabulationData tables(params);

    aggregator.tabulateNewManyBodyGrid(tables);

    auto& fe_fe_grid = tables.eam_grids_vec.back().aggregator_grids.getValueGrid({"Fe", "Fe"});
    for (auto cell: fe_fe_grid) {
        EXPECT_NEAR(cell.value, 1.0, 1e-9);
    }

    auto& fe_ni_grid = tables.eam_grids_vec.back().aggregator_grids.getValueGrid({"Fe", "Ni"});
    for (auto cell: fe_ni_grid) {
        EXPECT_NEAR(cell.value, 10.0, 1e-9);
    }
}

TEST(TestTwoBodySum, RealPairFunction) {
    Atoms atoms({{0, 0, 0}, {2, 0, 0}}, {Species("Fe"), Species("Ni")});
    auto nl = NeighbourLists(atoms, 5.0);

    auto aggregator = TwoBodySum<1>("Fe");
    aggregator.extend({"Fe", "Ni"}, PolycutoffPairFunction(4.0, 1.0, 2.0));

    auto aggregated_descriptors = aggregator.aggregate(nl);

    ASSERT_EQ(aggregated_descriptors.values.size(), 1);

    // PolycutoffPairFunction evaluation:
    // chi = (2.0 - 1.0) / (4.0 - 1.0) = 1/3
    // val = 2.0 * (1 - (1/27) * (6(1/9) - 15(1/3) + 10)) = 128/81
    Real expected_val = 128.0 / 81.0;
    EXPECT_NEAR(aggregated_descriptors.values[0][0], expected_val, 1e-9);

    // deriv = 2.0 * (1/3) * (1/9) * (-30(1/9) + 60(1/3) - 30) = -80/81
    Real expected_deriv = -80.0 / 81.0;

    // Force on atom 0: +direction * deriv = (1,0,0) * (-80/81)
    const auto& calc_force0 = aggregated_descriptors.force(0, 0)[0];
    EXPECT_NEAR(calc_force0.x, expected_deriv, 1e-9);
    EXPECT_NEAR(calc_force0.y, 0.0, 1e-9);
    EXPECT_NEAR(calc_force0.z, 0.0, 1e-9);

    // Virials: V += sep.virials * deriv
    // sep.virials.xx = -r_x*u_x = -2.0 * 1.0 = -2.0
    const auto& calc_virials0 = aggregated_descriptors.virials[0][0];
    EXPECT_NEAR(calc_virials0.xx, -2.0 * expected_deriv, 1e-9);
}

TEST(TestTwoBodySum, NonRotationallyInvariantThrowsOnExtend) {
    auto aggregator = TwoBodySum<1>("Fe");
    EXPECT_THROW(
        aggregator.extend({"Fe", "Ni"}, MockEamPairFunction(1.0, 0.5, 10.0, /*is_rot_inv=*/false)), std::runtime_error
    );
}
