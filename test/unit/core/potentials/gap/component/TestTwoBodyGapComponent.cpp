#include <gtest/gtest.h>
#include "jgap/core/ValuePtr.hpp"
#include "jgap/core/atomic/Atoms.hpp"
#include "jgap/core/atomic/neighbours/NeighbourLists.hpp"
#include "jgap/core/atomic/species/composition/Species2Sorted.hpp"
#include "jgap/core/cutoff/CosCutoff.hpp"
#include "jgap/core/kernels/Kernel.hpp"
#include "jgap/core/kernels/SquaredExpKernel.hpp"
#include "jgap/core/potentials/gap/component/TwoBodyGapComponent.hpp"
#include "jgap/core/tabulation/TabulationData.hpp"
#include "jgap/core/transform/nbody/2b/PairDistanceTransformation.hpp"
#include "jgap/core/transform/nbody/2b/TwoBodyTransformation.hpp"

using namespace jgap;

namespace {
    template<size_t Dim>
    class MockTwoBodyTransformation : public TwoBodyTransformation<Dim> {
    public:
        MockTwoBodyTransformation(bool is_rot_inv = true) : is_rot_inv(is_rot_inv) {}
        TwoBodyDescriptor<Dim> evaluateAndDifferentiate(const Cluster2& cluster) const override {
            TwoBodyDescriptor<Dim> res;
            res.value[0] = 1.0;
            const auto& dir = cluster.separation01.direction;
            for (size_t d = 0; d < Dim; ++d) {
                res.grad_r1[d] = 0.1_r * dir;
            }
            return res;
        }
        Cutoffs getCutoffs() const override { return Cutoffs{ { 2, 15.0 } }; }
        bool isRotationallyInvariant() const override { return is_rot_inv; }
        Descriptor<Dim> evaluate(const Cluster2& cluster) const override {
            return TwoBodyTransformation<Dim>::evaluate(cluster);
        }
        MockTwoBodyTransformation* clone() const override { return new MockTwoBodyTransformation(*this); }
    private:
        bool is_rot_inv = true;
    };

    template<size_t Dim>
    class MockKernel : public Kernel<Dim> {
    public:
        using KernelValueAndGradient = typename Kernel<Dim>::KernelValueAndGradient;
        Real value(const Descriptor<Dim>& q1, const Descriptor<Dim>& q2) const override {
            return Kernel<Dim>::value(q1, q2);
        }
        KernelValueAndGradient valueAndGradient(
            const Descriptor<Dim>& sparse_point, const Descriptor<Dim>& q
        ) const override {
            KernelValueAndGradient res;
            res.value = 2.0;
            res.gradient.fill(0.5);
            return res;
        }

        Kernel<Dim>* clone() const override { return new MockKernel<Dim>(); }
    };
}

TEST(TestTwoBodyGapComponent, MockTwoBody) {
    Atoms atoms(
        { { 0, 0, 0 }, { 10, 0, 0 }, { 0, 10, 0 }, { 10, 10, 0 } },
        { Species("Fe"), Species("Ni"), Species("Fe"), Species("Ni") }
    );
    auto nl = NeighbourLists(atoms, 11.0);

    auto component = TwoBodyGapComponent<1, MockKernel<1>>(
        Species2Sorted("Fe", "Ni"), MockTwoBodyTransformation<1>(), MockKernel<1>(), { { 1.0 } }
    );

    auto result = component.covariate(nl);
    ASSERT_TRUE(result.has_value());
    auto& quantities = result.value();

    // 2 Fe-Ni pairs. K=2.0. Energy = 2 * 2.0 = 4.0. Factor = 2.0 -> 8.0.
    EXPECT_NEAR(quantities.energy(0), 8.0, 1e-9);

    // Forces: dK/dr = (dq/dr) * (dK/dq) = 0.1 * 0.5 = 0.05. * 2 = 0.1.
    // F_i = dK/dr * direction
    Vector3 force0{ 0, 0, 0 }, force1{ 0, 0, 0 }, force2{ 0, 0, 0 }, force3{ 0, 0, 0 };
    Virials virials{};

    // Pair (0,1): Fe-Ni. direction = (1,0,0)
    force0 += Vector3{ 1, 0, 0 } * 0.1;
    force1 -= Vector3{ 1, 0, 0 } * 0.1;
    virials += Separation(atoms.getPositions()[0], atoms.getPositions()[1]).virials() * 0.1;

    // Pair (2,3): Fe-Ni. direction = (1,0,0)
    force2 += Vector3{ 1, 0, 0 } * 0.1;
    force3 -= Vector3{ 1, 0, 0 } * 0.1;
    virials += Separation(atoms.getPositions()[2], atoms.getPositions()[3]).virials() * 0.1;

    EXPECT_NEAR(quantities.force(0, 0).x, force0.x, 1e-9);
    EXPECT_NEAR(quantities.force(0, 1).x, force1.x, 1e-9);
    EXPECT_NEAR(quantities.force(0, 2).x, force2.x, 1e-9);
    EXPECT_NEAR(quantities.force(0, 3).x, force3.x, 1e-9);
    EXPECT_NEAR(quantities.virials(0).xx, virials.xx, 1e-9);
}

TEST(TestTwoBodyGapComponent, RealTwoBody) {
    Atoms atoms({ { 0, 0, 0 }, { 4, 0, 0 } }, { Species("Fe"), Species("Ni") });
    auto nl = NeighbourLists(atoms, 6.0);
    ValuePtr<TwoBodyTransformation<2>> trans = PairDistanceTransformation(CosCutoff(5.0, 2.0));
    auto kernel = SquaredExpKernel<1, 1>(1.0, std::array{ 1.0_r });
    Real expected_cutoff_val = 0.5;
    std::vector<Descriptor<2>> sparse_points = { { 4.0, expected_cutoff_val } };
    auto component = TwoBodyGapComponent(Species2Sorted("Fe", "Ni"), trans, kernel, sparse_points);

    auto result = component.covariate(nl);
    ASSERT_TRUE(result.has_value());
    auto& quantities = result.value();

    // (i,j) + (j,i) => x2
    EXPECT_NEAR(quantities.energy(0), 2.0 * expected_cutoff_val * expected_cutoff_val, 1e-9);

    // ~value drops slower as cutoff increases
    EXPECT_NEAR(quantities.force(0, 0).x, -0.78539816339744828, 1e-9);
    EXPECT_NEAR(quantities.force(0, 1).x, 0.78539816339744828, 1e-9);
    EXPECT_NEAR(quantities.virials(0).xx, 3.1415926535897931, 1e-9);
}

TEST(TestTwoBodyGapComponent, TabulationTwoBody) {
    ValuePtr<TwoBodyTransformation<2>> trans = PairDistanceTransformation(CosCutoff(5.0, 2.0));
    auto kernel = SquaredExpKernel<1, 1>(1.0, std::array{ 1.0_r });
    std::vector<Descriptor<2>> sparse_points = { { 4.0, 0.5 }, { 3.5, 0.5 } };
    std::vector<Real> coeffs = { 2.0, -1.0 };
    Species2Sorted species{ "Fe", "Ni" };
    auto component = TwoBodyGapComponent(species, trans, kernel, sparse_points);

    TabulationParams params;
    params.n_grid_2b = 10;
    params.max_cutoffs = component.getCutoffs();
    TabulationData tables(params);

    ASSERT_THROW(component.tabulate(tables), std::runtime_error);

    component.setCoefficients(coeffs);
    component.tabulate(tables);

    auto& grid = tables.two_body_grids.getValueGrid(species);

    for (size_t i = 0; i < params.n_grid_2b; i++) {
        Real r = grid.getCoord({ i })[0];

        Real K1 = std::exp(-0.5 * std::pow(r - 4.0, 2)) * CosCutoff(5.0, 2.0).evaluate(r) * 0.5;
        Real K2 = std::exp(-0.5 * std::pow(r - 3.5, 2)) * CosCutoff(5.0, 2.0).evaluate(r) * 0.5;

        Real expected_energy = 2.0 * (coeffs[0] * K1 + coeffs[1] * K2);

        EXPECT_NEAR(grid.data_flat[i], expected_energy, 1e-9);
    }
}

TEST(TestTwoBodyGapComponent, NonRotationallyInvariantThrowsOnTabulation) {
    ValuePtr<TwoBodyTransformation<1>> non_rot_trans = MockTwoBodyTransformation<1>(/*is_rot_inv=*/false);
    MockKernel<1> kernel;
    std::vector<Descriptor<1>> sparse_points = { { 1.0 } };
    Species2Sorted species{ "Fe", "Ni" };
    auto component = TwoBodyGapComponent(species, non_rot_trans, kernel, sparse_points);
    component.setCoefficients({ 1.0 });

    TabulationParams params;
    params.n_grid_2b = 10;
    params.max_cutoffs = component.getCutoffs();
    TabulationData tables(params);

    EXPECT_THROW(component.tabulate(tables), std::runtime_error);
}
