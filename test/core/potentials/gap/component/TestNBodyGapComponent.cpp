#include <gtest/gtest.h>
#include "core/potentials/gap/component/NBodyGapComponent.hpp"
#include "core/atomic/Atoms.hpp"
#include "core/atomic/neighbours/NeighbourList.hpp"
#include "core/transform/ClusterTransformation.hpp"
#include "core/kernels/Kernel.hpp"
#include "core/transform/2b/TwoBodyTransformation.hpp"
#include "core/transform/3b/Angle3bTransformation.hpp"
#include "core/kernels/SquaredExpKernel.hpp"
#include "core/cutoff/CosCutoff.hpp"

using namespace jgap;

namespace {
    template<size_t Dim, size_t ClusterSize>
    class MockClusterTransformation : public ClusterTransformation<Dim, ClusterSize> {
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
        Real symmetryFactor() const override {
            return ClusterSize == 3 ? 2.0 : 1.0;
        }

        std::unique_ptr<ClusterTransformation<Dim, ClusterSize>> clone() const override {
            return std::make_unique<MockClusterTransformation<Dim, ClusterSize>>();
        }
    };

    template<size_t Dim>
    class MockKernel : public Kernel<Dim> {
    public:
        using KernelValueAndGradient = Kernel<Dim>::KernelValueAndGradient;
        Real value(const std::array<Real, Dim>& q1, const std::array<Real, Dim>& q2) const override {
            return 2.0;
        }
        KernelValueAndGradient valueAndGradient(const std::array<Real, Dim>& sparse_point, const std::array<Real, Dim>& q) const override {
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

    auto component = NBodyGapComponent<1, 2, HasCentralAtom, MockKernel<1>>(
        SpeciesSet<2, HasCentralAtom>{"Fe", "Ni"},
        MockClusterTransformation<1, 2>(),
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

    auto component = NBodyGapComponent<1, 3, HasCentralAtom, MockKernel<1>>(
        SpeciesSet<3, HasCentralAtom>{"Fe", "Fe", "Ni"},
        MockClusterTransformation<1, 3>(),
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
    Atoms atoms({ {0,0,0}, {3,4,0} }, { Species("Fe"), Species("Ni") });
    auto nl = NeighbourList(atoms, 6.0);
    auto trans = TwoBodyTransformation(CosCutoff(5.0, 2.0));
    auto kernel = SquaredExpKernel<1, 1>(1.0, std::array{1.0});
    std::vector<Descriptor<2>> sparse_points = { {{{4.0, 0.5}}} };
    auto component = NBodyGapComponent<2, 2, Symmetric, SquaredExpKernel<1, 1>>(
        SpeciesSet<2, Symmetric>{"Fe", "Ni"}, trans, kernel, sparse_points
        );

    auto result = component.covariate(nl);
    ASSERT_TRUE(result.has_value());
    auto& quantities = result.value();

    // Pair (0,1) is at r=5.0, which is the cutoff. f_cut=0, so K=0.
    EXPECT_NEAR(quantities.energy(0), 0.0, 1e-9);
    EXPECT_NEAR(quantities.force(0, 0).norm(), 0.0, 1e-9);
    EXPECT_NEAR(quantities.force(0, 1).norm(), 0.0, 1e-9);
    EXPECT_NEAR(quantities.virials(0).xx, 0.0, 1e-9);
}

TEST(TestNBodyGapComponent, RealThreeBody) {
    Atoms atoms({ {0,0,0}, {3,0,0}, {0,4,0} }, { Species("Fe"), Species("Fe"), Species("Ni") });
    auto nl = NeighbourList(atoms, 6.0);
    auto trans = Angle3bTransformation(CosCutoff(5.0, 2.0));
    auto kernel = SquaredExpKernel<3, 1>(1.0, std::array{1.0, 1.0, 1.0});
    std::vector<Descriptor<4>> sparse_points = { {{{7.0, 1.0, 5.0, 0.25}}} };
    auto component = NBodyGapComponent<4, 3, HasCentralAtom, SquaredExpKernel<3, 1>>(
        SpeciesSet<3, HasCentralAtom>{"Fe", "Fe", "Ni"}, trans, kernel, sparse_points
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