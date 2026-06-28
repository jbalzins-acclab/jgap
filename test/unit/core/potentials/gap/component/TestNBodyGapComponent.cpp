#include <gtest/gtest.h>
#include "core/potentials/gap/component/NBodyGapComponent.hpp"
#include "core/atomic/Atoms.hpp"
#include "core/atomic/neighbours/NeighbourList.hpp"
#include "core/transform/NBodyTransformation.hpp"
#include "core/kernels/Kernel.hpp"
#include "core/transform/2b/TwoBodyTransformation.hpp"
#include "core/transform/3b/Angle3bTransformation.hpp"
#include "core/kernels/SquaredExpKernel.hpp"
#include "core/cutoff/CosCutoff.hpp"
#include "core/tabulation/TabulationData.hpp"

using namespace jgap;

namespace {
    template<size_t Dim, size_t ClusterSize>
    class MockTransformation : public NBodyTransformation<Dim, ClusterSize> {
    public:
        NBodyDescriptor<Dim, ClusterSize> evaluateAndDifferentiate(const Cluster<ClusterSize>& cluster)
            const override {
            NBodyDescriptor<Dim, ClusterSize> res;
            res.value[0] = 1.0;
            for (auto& deriv_dim : res.derivatives) {
                deriv_dim.fill(0.1);
            }
            return res;
        }
        Cutoffs getCutoffs() const override { return Cutoffs{{ClusterSize, 15.0}}; }
        Descriptor<Dim> evaluate(const Cluster<ClusterSize> &cluster) const override {
            throw std::runtime_error("Not implemented");
        }
        bool isSymmetricWrtNodeIndices() const override {
            return true;
        }
        Real symmetryFactor() const override {
            return ClusterSize == 3 ? 2.0 : 1.0;
        }

        MockTransformation* clone() const override {
            return new MockTransformation();
        }
    };

    template<size_t Dim>
    class MockKernel : public Kernel<Dim> {
    public:
        using KernelValueAndGradient = Kernel<Dim>::KernelValueAndGradient;
        Real value(const Descriptor<Dim>& q1, const Descriptor<Dim>& q2) const override {
            return 2.0;
        }
        KernelValueAndGradient valueAndGradient(const Descriptor<Dim>& sparse_point,
                                                const Descriptor<Dim>& q) const override {
            KernelValueAndGradient res;
            res.value = 2.0;
            res.gradient.fill(0.5);
            return res;
        }
    };
}

TEST(TestNBodyGapComponent, MockTwoBody) {
    Atoms atoms({ {0,0,0}, {10,0,0}, {0,10,0}, {10,10,0} },
        { Species("Fe"), Species("Ni"), Species("Fe"), Species("Ni") });
    auto nl = NeighbourList(atoms, 11.0);

    auto component = NBodyGapComponent<1, 2, NodeSymmetric, MockKernel<1>>(
        SpeciesSet<2, NodeSymmetric>{"Fe", "Ni"},
        MockTransformation<1, 2>(),
        MockKernel<1>(),
        { {{{1.0}}} }
    );

    auto result = component.covariate(nl);
    ASSERT_TRUE(result.has_value());
    auto& quantities = result.value();

    // Energy: 2 Fe-Ni pairs * (K=2.0) = 4.0
    EXPECT_NEAR(quantities.energy(0), 4.0, 1e-9);

    // Forces: dK/dr = (dq/dr) * (dK/dq) = 0.1 * 0.5 = 0.05
    // F_i = dK/dr * direction
    Vector3 force0{0,0,0}, force1{0,0,0}, force2{0,0,0}, force3{0,0,0};
    Virials virials{};

    // Pair (0,1): Fe-Ni. direction = (1,0,0)
    force0 += Vector3{1,0,0} * 0.05;
    force1 -= Vector3{1,0,0} * 0.05;
    virials += Separation(atoms.lookupPositions()[0], atoms.lookupPositions()[1]).derivatives.virials * 0.05;

    // Pair (2,3): Fe-Ni. direction = (1,0,0)
    force2 += Vector3{1,0,0} * 0.05;
    force3 -= Vector3{1,0,0} * 0.05;
    virials += Separation(atoms.lookupPositions()[2], atoms.lookupPositions()[3]).derivatives.virials * 0.05;

    EXPECT_NEAR(quantities.force(0, 0).x, force0.x, 1e-9);
    EXPECT_NEAR(quantities.force(0, 1).x, force1.x, 1e-9);
    EXPECT_NEAR(quantities.force(0, 2).x, force2.x, 1e-9);
    EXPECT_NEAR(quantities.force(0, 3).x, force3.x, 1e-9);
    EXPECT_NEAR(quantities.virials(0).xx, virials.xx, 1e-9);
}

TEST(TestNBodyGapComponent, MockThreeBody) {
    Atoms atoms({ {0,0,0}, {10,0,0}, {0,10,0}, {10,10,0} }, { Species("Fe"), Species("Ni"), Species("Fe"), Species("Ni") });
    auto nl = NeighbourList(atoms, 11.0);

    auto component = NBodyGapComponent<1, 3, NodeSymmetric, MockKernel<1>>(
        SpeciesSet<3, NodeSymmetric>{"Fe", "Fe", "Ni"},
        MockTransformation<1, 3>(),
        MockKernel<1>(),
        { {{{1.0}}} }
    );

    auto result = component.covariate(nl);
    ASSERT_TRUE(result.has_value());
    auto& quantities = result.value();

    // Energy: 2 triplets * (K=2.0 * sym=2.0) = 8.0
    EXPECT_NEAR(quantities.energy(0), 8.0, 1e-9);

    // Forces: dK/dr = (dq/dr * dK/dq) * sym = (0.1 * 0.5) * 2.0 = 0.1
    Vector3 force0{0,0,0}, force1{0,0,0}, force2{0,0,0}, force3{0,0,0};
    Virials virials{};

    // Triplet (0,2,1)
    Separation sep02(atoms.lookupPositions()[0], atoms.lookupPositions()[2]);
    force0 += sep02.derivatives.direction * 0.1;
    force2 -= sep02.derivatives.direction * 0.1;
    virials += sep02.derivatives.virials * 0.1;
    Separation sep01(atoms.lookupPositions()[0], atoms.lookupPositions()[1]);
    force0 += sep01.derivatives.direction * 0.1;
    force1 -= sep01.derivatives.direction * 0.1;
    virials += sep01.derivatives.virials * 0.1;
    Separation sep21(atoms.lookupPositions()[2], atoms.lookupPositions()[1]);
    force2 += sep21.derivatives.direction * 0.1;
    force1 -= sep21.derivatives.direction * 0.1;
    virials += sep21.derivatives.virials * 0.1;

    // Triplet (2,0,3)
    Separation sep20(atoms.lookupPositions()[2], atoms.lookupPositions()[0]);
    force2 += sep20.derivatives.direction * 0.1;
    force0 -= sep20.derivatives.direction * 0.1;
    virials += sep20.derivatives.virials * 0.1;
    Separation sep23(atoms.lookupPositions()[2], atoms.lookupPositions()[3]);
    force2 += sep23.derivatives.direction * 0.1;
    force3 -= sep23.derivatives.direction * 0.1;
    virials += sep23.derivatives.virials * 0.1;
    Separation sep03(atoms.lookupPositions()[0], atoms.lookupPositions()[3]);
    force0 += sep03.derivatives.direction * 0.1;
    force3 -= sep03.derivatives.direction * 0.1;
    virials += sep03.derivatives.virials * 0.1;

    EXPECT_NEAR(quantities.force(0, 0).x, force0.x, 1e-9);
    EXPECT_NEAR(quantities.force(0, 1).x, force1.x, 1e-9);
    EXPECT_NEAR(quantities.force(0, 2).x, force2.x, 1e-9);
    EXPECT_NEAR(quantities.force(0, 3).x, force3.x, 1e-9);
    EXPECT_NEAR(quantities.virials(0).xx, virials.xx, 1e-9);
}

TEST(TestNBodyGapComponent, RealTwoBody) {
    Atoms atoms({ {0,0,0}, {4,0,0} }, { Species("Fe"), Species("Ni") });
    auto nl = NeighbourList(atoms, 6.0);
    auto trans = TwoBodyTransformation(CosCutoff(5.0, 2.0));
    auto kernel = SquaredExpKernel<1, 1>(1.0, std::array{1.0});
    std::vector<Descriptor<2>> sparse_points = { {{{4.0, 1.0}}} };
    auto component = NBodyGapComponent<2, 2, FullSymmetry, SquaredExpKernel<1, 1>>(
        SpeciesSet<2, FullSymmetry>{"Fe", "Ni"}, trans, kernel, sparse_points
    );

    auto result = component.covariate(nl);
    ASSERT_TRUE(result.has_value());
    auto& quantities = result.value();

    // Pair (0,1) is at r=4.0. Cutoff is CosCutoff(5.0, 2.0). f_cut = 0.5 * (cos(pi * (4-3)/2) + 1) = 0.5
    // q = {4.0, 0.5}. qs = {4.0, 1.0}.
    // K = prefactor * exp(-0.5 * (4.0-4.0)^2 / 1.0) * (0.5 * 1.0) = 0.5
    // With sym_factor = 2.0, Energy = 0.5 * 2.0 = 1.0
    EXPECT_NEAR(quantities.energy(0), 1.0, 1e-9);

    // dK/dq[0] = 0. dK/dq[1] = 1.0. With sym_factor = 2.0, grad_q = {0.0, 2.0}
    // dq/dr = {1.0, -pi/4}. dK/drij = dq/dr * grad_q = -pi/2
    // direction_01 = {1.0, 0, 0}. force0 = {-pi/2, 0, 0}, force1 = {pi/2, 0, 0}
    EXPECT_NEAR(quantities.force(0, 0).x, -M_PI / 2.0, 1e-9);
    EXPECT_NEAR(quantities.force(0, 1).x, M_PI / 2.0, 1e-9);
    EXPECT_NEAR(quantities.force(0, 0).y, 0.0, 1e-9);
    EXPECT_NEAR(quantities.force(0, 1).y, 0.0, 1e-9);

    // virials = separation x force. separation = {4,0,0}. virials.xx = 4.0 * (pi/2) = 2*pi
    EXPECT_NEAR(quantities.virials(0).xx, 2.0 * M_PI, 1e-9);
}

TEST(TestNBodyGapComponent, RealThreeBody) {
    Atoms atoms({ {0,0,0}, {3,0,0}, {0,4,0} }, { Species("Fe"), Species("Fe"), Species("Ni") });
    auto nl = NeighbourList(atoms, 6.0);
    auto trans = Angle3bTransformation(CosCutoff(5.0, 2.0));
    auto kernel = SquaredExpKernel<3, 1>(1.0, std::array{1.0, 1.0, 1.0});
    std::vector<Descriptor<4>> sparse_points = { {{{7.0, 1.0, 5.0, 0.25}}} };
    auto component = NBodyGapComponent<4, 3, NodeSymmetric, SquaredExpKernel<3, 1>>(
        SpeciesSet<3, NodeSymmetric>{"Fe", "Fe", "Ni"}, trans, kernel, sparse_points
        );

    auto result = component.covariate(nl);
    ASSERT_TRUE(result.has_value());
    auto& quantities = result.value();

    // Only triplet (0,1,2) contributes, K = 0.125. SymFactor = 2.0. Total E = 0.25.
    EXPECT_NEAR(quantities.energy(0), 0.25, 1e-9);

    // dK/dr02 = -pi/8. F_0 = dK/dr02 * direction_02. F_2 = -dK/dr02 * direction_02
    Vector3 force0{0,0,0}, force1{0,0,0}, force2{0,0,0};
    Virials virials{};

    Separation sep02(atoms.lookupPositions()[0], atoms.lookupPositions()[2]);
    Real dK_dr02 = -M_PI / 8.0;
    force0 += sep02.derivatives.direction * dK_dr02;
    force2 -= sep02.derivatives.direction * dK_dr02;
    virials += sep02.derivatives.virials * dK_dr02;

    EXPECT_NEAR(quantities.force(0, 0).norm(), force0.norm(), 1e-9);
    EXPECT_NEAR(quantities.force(0, 1).norm(), force1.norm(), 1e-9);
    EXPECT_NEAR(quantities.force(0, 2).norm(), force2.norm(), 1e-9);
    EXPECT_NEAR(quantities.virials(0).xx, virials.xx, 1e-9);
}

TEST(TestNBodyGapComponent, TabulationTwoBody) {
    auto trans = TwoBodyTransformation(CosCutoff(5.0, 2.0));
    auto kernel = SquaredExpKernel<1, 1>(1.0, std::array{1.0});
    std::vector<Descriptor<2>> sparse_points = { {{{4.0, 1.0}}}, {{{3.5, 0.5}}} };
    std::vector<Real> coeffs = {2.0, -1.0};
    SpeciesSet<2, FullSymmetry> species{"Fe", "Ni"};
    auto component = NBodyGapComponent<2, 2, FullSymmetry, SquaredExpKernel<1, 1>>(
        species, trans, kernel, sparse_points
    );

    TabulationParams params;
    params.n_grid_2b = 10;
    params.max_cutoffs = component.getCutoffs();
    TabulationData tables(params);

    ASSERT_THROW(component.tabulate(tables), std::runtime_error);

    component.setCoefficients(coeffs);
    component.tabulate(tables);

    auto& grid = tables.two_body_grids.getValueGrid(species);
    auto cutoff = CosCutoff(5.0, 2.0);

    for (size_t i = 0; i < params.n_grid_2b; i++) {
        Real r = grid.getCoord({i})[0];
        Real f_cut = cutoff.evaluate(r);

        Real K1 = std::exp(-0.5 * std::pow(r - 4.0, 2)) * (f_cut * 1.0);
        Real K2 = std::exp(-0.5 * std::pow(r - 3.5, 2)) * (f_cut * 0.5);

        Real expected_energy = 2.0 * (coeffs[0] * K1 + coeffs[1] * K2);

        EXPECT_NEAR(grid.data_flat[i], expected_energy, 1e-9);
    }
}

TEST(TestNBodyGapComponent, TabulationThreeBody) {
    auto trans = Angle3bTransformation(CosCutoff(5.0, 2.0));
    auto kernel = SquaredExpKernel<3, 1>(1.0, std::array{1.0, 1.0, 1.0});
    std::vector<Descriptor<4>> sparse_points = { {{{7.0, 1.0, 5.0, 0.25}}} };
    std::vector<Real> coeffs = {0.5};
    SpeciesSet<3, NodeSymmetric> species{"Fe", "Fe", "Ni"};
    auto component = NBodyGapComponent<4, 3, NodeSymmetric, SquaredExpKernel<3, 1>>(
        species, trans, kernel, sparse_points
    );

    TabulationParams params;
    params.n_grid_3b = {5, 5, 5};
    params.max_cutoffs = component.getCutoffs();
    TabulationData tables(params);

    ASSERT_THROW(component.tabulate(tables), std::runtime_error);

    component.setCoefficients(coeffs);
    component.tabulate(tables);

    auto& grid = tables.three_body_grids.getValueGrid(species);
    auto cutoff = CosCutoff(5.0, 2.0);

    for (size_t i = 0; i < grid.data_flat.size(); ++i) {
        auto grid_pos = grid.getCoord(grid.getIndices(i));
        auto triplet = TabulationData::gridPosAsCluster<3>(grid_pos);

        Real r01 = triplet.separationBetween(0, 1);
        Real r02 = triplet.separationBetween(0, 2);
        Real r12 = triplet.separationBetween(1, 2);

        Real f_cut_01 = cutoff.evaluate(r01);
        Real f_cut_02 = cutoff.evaluate(r02);

        auto q = std::array{
            r01 + r02,
            (r01 - r02) * (r01 - r02),
            r12,
            f_cut_01 * f_cut_02
        };

        auto qs = sparse_points[0].value;
        Real K = std::exp(-0.5 * (
            std::pow(q[0] - qs[0], 2) / 1.0 +
            std::pow(q[1] - qs[1], 2) / 1.0 +
            std::pow(q[2] - qs[2], 2) / 1.0
        )) * (q[3] * qs[3]);

        EXPECT_NEAR(grid.data_flat[i], 2.0 * coeffs[0] * K, 1e-9);
    }
}