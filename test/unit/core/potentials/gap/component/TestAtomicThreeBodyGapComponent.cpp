#include <gtest/gtest.h>
#include "jgap/core/ValuePtr.hpp"
#include "jgap/core/atomic/Atoms.hpp"
#include "jgap/core/atomic/neighbours/NeighbourLists.hpp"
#include "jgap/core/atomic/species/composition/Species3AtomicSorted.hpp"
#include "jgap/core/cutoff/CosCutoff.hpp"
#include "jgap/core/kernels/Kernel.hpp"
#include "jgap/core/kernels/SquaredExpKernel.hpp"
#include "jgap/core/potentials/gap/component/AtomicThreeBodyGapComponent.hpp"
#include "jgap/core/tabulation/TabulationData.hpp"
#include "jgap/core/transform/nbody/3b/Angle3bTransformation.hpp"
#include "jgap/core/transform/nbody/3b/ThreeBodyTransformation.hpp"

using namespace jgap;

namespace {
    template<size_t Dim>
    class MockThreeBodyTransformation : public ThreeBodyTransformation<Dim> {
    public:
        ThreeBodyDescriptor<Dim> evaluateAndDifferentiate(const Cluster3& cluster) const override {
            ThreeBodyDescriptor<Dim> res;
            res.value[0] = 1.0;
            res.derivatives[0].fill(0.1);
            res.derivatives[1].fill(0.1);
            res.derivatives[2].fill(0.1);
            return res;
        }
        Cutoffs getCutoffs() const override { return Cutoffs{{3, 15.0}}; }
        Descriptor<Dim> evaluate(const Cluster3& cluster) const override {
            Descriptor<Dim> res;
            res[0] = 1.0;
            return res;
        }
        bool isSwapInvariant(size_t idx1, size_t idx2) const override {
            return (idx1 == 1 && idx2 == 2) || (idx1 == 2 && idx2 == 1);
        }
        MockThreeBodyTransformation* clone() const override { return new MockThreeBodyTransformation(); }
    };

    template<size_t Dim>
    class MockKernel : public Kernel<Dim> {
    public:
        using KernelValueAndGradient = typename Kernel<Dim>::KernelValueAndGradient;
        Real value(const Descriptor<Dim>& q1, const Descriptor<Dim>& q2) const override { return 2.0; }
        KernelValueAndGradient valueAndGradient(const Descriptor<Dim>& sparse_point,
                                                const Descriptor<Dim>& q) const override {
            KernelValueAndGradient res;
            res.value = 2.0;
            res.gradient.fill(0.5);
            return res;
        }
    };
}

TEST(TestAtomicThreeBodyGapComponent, MockThreeBody) {
    Atoms atoms({{0, 0, 0}, {10, 0, 0}, {0, 10, 0}, {10, 10, 0}},
                {Species("Fe"), Species("Ni"), Species("Fe"), Species("Ni")});
    auto nl = NeighbourLists(atoms, 11.0);

    auto component = AtomicThreeBodyGapComponent<1, MockKernel<1>>(
        Species3AtomicSorted("Fe", "Fe", "Ni"), MockThreeBodyTransformation<1>(), MockKernel<1>(), {{1.0}});

    auto result = component.covariate(nl);
    ASSERT_TRUE(result.has_value());
    auto& quantities = result.value();

    // 2 triplets centered at root Fe atoms: (0,2,1) and (2,0,3).
    // Each triplet has K = 2.0.
    // iteration_reduction_factor = 2.0 because isSwapInvariant(1,2) == true.
    // Energy: 2 triplets * K=2.0 * reduction_factor=2.0 = 8.0
    EXPECT_NEAR(quantities.energy(0), 8.0, 1e-9);

    // Forces: dK/dr = (dq/dr * dK/dq) = 0.1 * 0.5 = 0.05. * 2.0 = 0.1
    Vector3 force0{0, 0, 0}, force1{0, 0, 0}, force2{0, 0, 0}, force3{0, 0, 0};
    Virials virials{};

    // Triplet (0,2,1)
    Separation sep02(atoms.getPositions()[0], atoms.getPositions()[2]);
    force0 += sep02.derivatives.direction * 0.1;
    force2 -= sep02.derivatives.direction * 0.1;
    virials += sep02.derivatives.virials * 0.1;

    Separation sep01(atoms.getPositions()[0], atoms.getPositions()[1]);
    force0 += sep01.derivatives.direction * 0.1;
    force1 -= sep01.derivatives.direction * 0.1;
    virials += sep01.derivatives.virials * 0.1;

    Separation sep21(atoms.getPositions()[2], atoms.getPositions()[1]);
    force2 += sep21.derivatives.direction * 0.1;
    force1 -= sep21.derivatives.direction * 0.1;
    virials += sep21.derivatives.virials * 0.1;

    // Triplet (2,0,3)
    Separation sep20(atoms.getPositions()[2], atoms.getPositions()[0]);
    force2 += sep20.derivatives.direction * 0.1;
    force0 -= sep20.derivatives.direction * 0.1;
    virials += sep20.derivatives.virials * 0.1;

    Separation sep23(atoms.getPositions()[2], atoms.getPositions()[3]);
    force2 += sep23.derivatives.direction * 0.1;
    force3 -= sep23.derivatives.direction * 0.1;
    virials += sep23.derivatives.virials * 0.1;

    Separation sep03(atoms.getPositions()[0], atoms.getPositions()[3]);
    force0 += sep03.derivatives.direction * 0.1;
    force3 -= sep03.derivatives.direction * 0.1;
    virials += sep03.derivatives.virials * 0.1;

    EXPECT_NEAR(quantities.force(0, 0).x, force0.x, 1e-9);
    EXPECT_NEAR(quantities.force(0, 1).x, force1.x, 1e-9);
    EXPECT_NEAR(quantities.force(0, 2).x, force2.x, 1e-9);
    EXPECT_NEAR(quantities.force(0, 3).x, force3.x, 1e-9);
    EXPECT_NEAR(quantities.virials(0).xx, virials.xx, 1e-9);
}

TEST(TestAtomicThreeBodyGapComponent, RealThreeBody) {
    Atoms atoms({{0, 0, 0}, {3, 0, 0}, {0, 4, 0}}, {Species("Fe"), Species("Fe"), Species("Ni")});
    auto nl = NeighbourLists(atoms, 6.0);
    ValuePtr<ThreeBodyTransformation<4>> trans = Angle3bTransformation(CosCutoff(5.0, 2.0));
    auto kernel = SquaredExpKernel<3, 1>(1.0, std::array{1.0, 1.0, 1.0});
    std::vector<Descriptor<4>> sparse_points = {{7.0, 1.0, 5.0, 0.25}};
    auto component = AtomicThreeBodyGapComponent(Species3AtomicSorted("Fe", "Fe", "Ni"), trans, kernel, sparse_points);

    auto result = component.covariate(nl);
    ASSERT_TRUE(result.has_value());
    auto& quantities = result.value();

    // Triplet (0,1,2) centered at Fe(0) with Fe(1) and Ni(2) as neighbors.
    // r01 = 3.0, r02 = 4.0, r12 = 5.0.
    // Cutoffs: f_cut(3.0) = 1.0, f_cut(4.0) = 0.5.
    // q = {7.0, 1.0, 5.0, 0.5}.
    // qs = {7.0, 1.0, 5.0, 0.25}.
    // K_raw = 0.5 * 0.25 = 0.125.
    // Angle3bTransformation is swap invariant (1,2), so iteration_reduction_factor = 2.0.
    // Total E = 0.25.
    EXPECT_NEAR(quantities.energy(0), 0.25, 1e-9);
}

TEST(TestAtomicThreeBodyGapComponent, TabulationThreeBody) {
    ValuePtr<ThreeBodyTransformation<4>> trans = Angle3bTransformation(CosCutoff(5.0, 2.0));
    auto kernel = SquaredExpKernel<3, 1>(1.0, std::array{1.0, 1.0, 1.0});
    std::vector<Descriptor<4>> sparse_points = {{7.0, 1.0, 5.0, 0.25}};
    std::vector<Real> coeffs = {0.5};
    Species3AtomicSorted species{"Fe", "Fe", "Ni"};
    auto component = AtomicThreeBodyGapComponent(species, trans, kernel, sparse_points);

    TabulationParams params;
    params.n_grid_3b = {5, 5, 5};
    params.max_cutoffs = component.getCutoffs();
    TabulationData tables(params);

    // Unset coefficients throw
    ASSERT_THROW(component.tabulate(tables), std::runtime_error);

    component.setCoefficients(coeffs);
    component.tabulate(tables);

    auto& grid = tables.three_body_grids.getValueGrid(species);
    EXPECT_GT(grid.data_flat.size(), 0);

    for (size_t i = 0; i < grid.data_flat.size(); ++i) {
        auto grid_pos = grid.getCoord(grid.getIndices(i));
        auto triplet = TabulationData::gridPosAsCluster3(grid_pos);

        auto transformed = trans->evaluate(triplet);
        Real K = kernel.value(sparse_points[0], transformed);

        EXPECT_NEAR(grid.data_flat[i], 2.0 * coeffs[0] * K, 1e-9);
    }
}
