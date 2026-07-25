#include <gtest/gtest.h>

#include "jgap/core/atomic/Atoms.hpp"
#include "jgap/core/atomic/io/XYZData.hpp"
#include "jgap/core/atomic/neighbours/NeighbourLists.hpp"
#include "jgap/core/atomic/species/composition/Species2Atomic.hpp"
#include "jgap/core/atomic/species/composition/Species2Sorted.hpp"
#include "jgap/core/atomic/species/composition/Species3AtomicSorted.hpp"
#include "jgap/core/cutoff/CosCutoff.hpp"
#include "jgap/core/fit/gap/regularization/SimpleRegularizationRules.hpp"
#include "jgap/core/kernels/SquaredExpKernel.hpp"
#include "jgap/core/potentials/gap/GapPotential.hpp"
#include "jgap/core/potentials/gap/component/AtomicThreeBodyGapComponent.hpp"
#include "jgap/core/potentials/gap/component/ManyBodyGapComponent.hpp"
#include "jgap/core/potentials/gap/component/TwoBodyGapComponent.hpp"
#include "jgap/core/transform/manybody/TwoBodySum.hpp"
#include "jgap/core/transform/nbody/2b/PairDistanceTransformation.hpp"
#include "jgap/core/transform/nbody/2b/eam/PolycutoffPairFunction.hpp"
#include "jgap/core/transform/nbody/3b/Angle3bTransformation.hpp"
#include "jgap/ext/fit/gap/QRGapFit.hpp"

using namespace jgap;

Atoms setupEquilateralTriangle() {
    return Atoms({{0.0, 0.0, 0.0}, {3.0, 0.0, 0.0}, {1.5, 2.598, 0.0}}, {Species("Fe"), Species("Fe"), Species("Fe")},
                 Lattice{{100.0, 0.0, 0.0}, {0.0, 100.0, 0.0}, {0.0, 0.0, 100.0}});
}

TEST(TestQRGapFit, twoBodyEquilateralTriangleAtEquilibriumQuipCompatibility) {
    auto equilateralTriangle = setupEquilateralTriangle();
    equilateralTriangle.setEnergy(1.0);
    equilateralTriangle.setForces({{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}});

    auto trans = PairDistanceTransformation(CosCutoff(10.0, 0.7));
    auto kernel = SquaredExpKernel<1, 1>(1.0, std::array<Real, 1>{1.0});
    std::vector<Descriptor<2>> sparse_points = {{3.0, 1.0}};
    auto component =
        TwoBodyGapComponent<2, SquaredExpKernel<1, 1>>(Species2Sorted("Fe", "Fe"), trans, kernel, sparse_points);

    auto potential = GapPotential({component});
    auto regularization_rules = SimpleRegularizationRules(3.0, 10.0, 5.0, 5.0);

    QRGapFit fitter;
    fitter.fit(potential, {equilateralTriangle}, regularization_rules);

    const auto& fitted_coeffs = potential.getComponents()[0]->getCoefficients();
    double c = fitted_coeffs[0];
    auto totalVariance = 6.0 / c - 36.0;
    ASSERT_NEAR(totalVariance, 27.0, 1e-6);
}

GapPotential create2bPotential(const std::vector<Real>& sparsePts) {
    auto trans = PairDistanceTransformation(CosCutoff(10.0, 0.7));
    auto kernel = SquaredExpKernel<1, 1>(1.0, std::array<Real, 1>{1.0});

    std::vector<Descriptor<2>> sparse_points;
    CosCutoff ref_cutoff(10.0, 0.7);
    for (Real r: sparsePts) {
        sparse_points.push_back({r, ref_cutoff.evaluate(r)});
    }

    auto component = TwoBodyGapComponent<2, SquaredExpKernel<1, 1>>(Species2Sorted("Fe", "Fe"), std::move(trans),
                                                                    kernel, sparse_points);

    return GapPotential({component});
}

Atoms setupTwoAtoms() {
    return Atoms({{0.0, 0.0, 0.0}, {3.0, 0.0, 0.0}}, {Species("Fe"), Species("Fe")},
                 Lattice{{100.0, 0.0, 0.0}, {0.0, 100.0, 0.0}, {0.0, 0.0, 100.0}});
}

TEST(TestQRGapFit, twoAtomsWithForceQuipCompatibility1) {
    auto twoAtoms = setupTwoAtoms();
    twoAtoms.setEnergy(1.0);
    twoAtoms.setForces({{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}});
    auto potential = create2bPotential({2.0, 4.0});
    auto rules = SimpleRegularizationRules(1.0, 1.0, 5.0, 5.0);
    QRGapFit fitter;
    fitter.fit(potential, {twoAtoms}, rules);
    const auto& coeffs = potential.getComponents()[0]->getCoefficients();
    ASSERT_NEAR(coeffs[0], 0.23266775926220731, 1e-6);
    ASSERT_NEAR(coeffs[1], 0.23266775926220731, 1e-6);
}

TEST(TestQRGapFit, twoAtomsWithForceQuipCompatibility2) {
    auto twoAtoms = setupTwoAtoms();
    twoAtoms.setEnergy(0.0);
    twoAtoms.setForces({{1.0, 0.0, 0.0}, {0.0, 0.0, 0.0}});
    auto potential = create2bPotential({2.0});
    auto rules = SimpleRegularizationRules(10000000.0, 1.0, 5.0, 5.0);
    QRGapFit fitter;
    fitter.fit(potential, {twoAtoms}, rules);
    const auto& coeffs = potential.getComponents()[0]->getCoefficients();
    ASSERT_NEAR(coeffs[0], -0.30764655994411466, 1e-6);
}

TEST(TestQRGapFit, twoAtomsWithForceQuipCompatibility3) {
    Atoms twoAtoms({{0.0, 0.0, 0.0}, {1.5, 2.598, 0.0}}, {Species("Fe"), Species("Fe")},
                   Lattice{{100.0, 0.0, 0.0}, {0.0, 100.0, 0.0}, {0.0, 0.0, 100.0}});
    twoAtoms.setEnergy(0.0);
    twoAtoms.setForces({{1.0, 0.0, 0.0}, {0.0, 0.0, 0.0}});
    auto potential = create2bPotential({2.0});
    auto rules = SimpleRegularizationRules(10000000.0, 1.0, 5.0, 5.0);
    QRGapFit fitter;
    fitter.fit(potential, {twoAtoms}, rules);
    const auto& coeffs = potential.getComponents()[0]->getCoefficients();
    ASSERT_NEAR(coeffs[0], -0.15382666452610533, 1e-6);
}

TEST(TestQRGapFit, twoAtomsWithForceQuipCompatibility4) {
    Atoms twoAtoms({{0.0, 0.0, 0.0}, {1.5, 2.598, 0.0}}, {Species("Fe"), Species("Fe")},
                   Lattice{{100.0, 0.0, 0.0}, {0.0, 100.0, 0.0}, {0.0, 0.0, 100.0}});
    twoAtoms.setEnergy(0.0);
    twoAtoms.setForces({{1.0, 1.0, 1.0}, {-1.0, -1.0, -1.0}});
    auto potential = create2bPotential({2.0});
    auto rules = SimpleRegularizationRules(10000000.0, 2.0, 5.0, 5.0);
    QRGapFit fitter;
    fitter.fit(potential, {twoAtoms}, rules);
    const auto& coeffs = potential.getComponents()[0]->getCoefficients();
    ASSERT_NEAR(coeffs[0], -0.47733536744854282, 1e-6);
}

TEST(TestQRGapFit, twoAtomsWithForceQuipCompatibility5) {
    Atoms twoAtoms({{0.0, 0.0, 0.0}, {1.5, 2.598, 0.0}}, {Species("Fe"), Species("Fe")},
                   Lattice{{100.0, 0.0, 0.0}, {0.0, 100.0, 0.0}, {0.0, 0.0, 100.0}});
    twoAtoms.setEnergy(0.0);
    twoAtoms.setForces({{1.0, 1.0, 1.0}, {-1.0, -1.0, -1.0}});
    auto potential = create2bPotential({2.0, 2.5, 4.0});
    auto rules = SimpleRegularizationRules(1.0, 2.0, 5.0, 5.0);
    QRGapFit fitter;
    fitter.fit(potential, {twoAtoms}, rules);
    const auto& coeffs = potential.getComponents()[0]->getCoefficients();
    ASSERT_NEAR(coeffs[0], -0.24567921601436607, 1e-6);
    ASSERT_NEAR(coeffs[1], -0.10998748761736651, 1e-6);
    ASSERT_NEAR(coeffs[2], 0.38697171383300527, 1e-6);
}

GapPotential createEamPotential(Real theta, Real delta, Real r_min, Real cutoff,
                                const std::vector<Real>& sparse_pts) {
    auto trans = PolycutoffPairFunction(cutoff, r_min, 1.0);

    auto eam_aggregator = TwoBodySum<1>(Species("Fe"));
    eam_aggregator.extend(Species2Atomic("Fe", "Fe"), trans);

    auto kernel = SquaredExpKernel<1, 0>(delta, std::array<Real, 1>{theta});

    std::vector<Descriptor<1>> sparse_points;
    for (Real r: sparse_pts) {
        sparse_points.push_back({r});
    }

    auto component = ManyBodyGapComponent<1, SquaredExpKernel<1, 0>>(eam_aggregator, kernel, sparse_points);

    return GapPotential({component});
}

TEST(TestQRGapFit, twoAtomsEamQuipCompatibility) {
    auto twoAtoms = setupTwoAtoms();
    twoAtoms.setEnergy(1.0);
    twoAtoms.setForces({{1.0, 0.0, 0.0}, {-1.0, 0.0, 0.0}});
    auto potential = createEamPotential(1.0, 1.0, 0.0, 5.0, {1.0});
    auto rules = SimpleRegularizationRules(1.0, 1.0, 1.0, 1.0);
    QRGapFit fitter;
    fitter.fit(potential, {twoAtoms}, rules);
    const auto& coeffs = potential.getComponents()[0]->getCoefficients();
    ASSERT_NEAR(coeffs[0], 0.017637586136984833, 1e-6);
}

TEST(TestQRGapFit, eamQuipCompatibilityRealBox) {
    auto box = Atoms::readAtoms("test/resources/xyz-samples/fe-only.xyz")[15];
    box.eraseVirials();
    auto potential = createEamPotential(3.0, 2.0, 0.0, 5.0, {1.0, 2.5, 4.0});
    auto rules = SimpleRegularizationRules(0.001, 0.05, 1.0, 1.0);
    QRGapFit fitter;
    fitter.fit(potential, {box}, rules);
    const auto& coeffs = potential.getComponents()[0]->getCoefficients();

    /*
        <sparseX i="1" alpha="-5.0619173053329485" sparseCutoff="1.0000000000000000"/>
        <sparseX i="2" alpha="-4.4169841885357615" sparseCutoff="1.0000000000000000"/>
        <sparseX i="3" alpha="-.16117464517539212" sparseCutoff="1.0000000000000000"/>
    */

    ASSERT_NEAR(coeffs[0], -5.0619173053329485, 1e-6);
    ASSERT_NEAR(coeffs[1], -4.4169841885357615, 1e-6);
    ASSERT_NEAR(coeffs[2], -0.16117464517539212, 1e-6);
}

GapPotential create3bPotential(Real theta, Real delta, Real cutoff_transition_width, Real cutoff,
                               const std::vector<Vector3>& sparsePts) {
    auto trans = Angle3bTransformation(CosCutoff(cutoff, cutoff_transition_width));
    auto kernel = SquaredExpKernel<3, 1>(delta, std::array<Real, 3>{theta, theta, theta});
    std::vector<Descriptor<4>> sparse_points;
    for (const auto& q: sparsePts) {
        sparse_points.push_back({q.x, q.y, q.z, 1.0});
    }
    auto component = AtomicThreeBodyGapComponent<4, SquaredExpKernel<3, 1>>(
        Species3AtomicSorted("Fe", "Fe", "Fe"), std::move(trans), std::move(kernel), sparse_points);

    return GapPotential({component});
}

TEST(TestQRGapFit, equilateralTriangle3bQuipCompatibility) {
    auto equilateralTriangle = setupEquilateralTriangle();
    equilateralTriangle.setEnergy(1.0);
    equilateralTriangle.setForces({{1.0, 0.0, 0.0}, {0.0, -1.0, 0.0}, {0.5, 0.5, 0.0}});
    auto potential = create3bPotential(1.0, 1.0, 0.6, 10.0, {{6.0, 0.0, 3.0}});
    auto rules = SimpleRegularizationRules(1.0, 1.0, 1.0, 1.0);
    QRGapFit fitter;
    fitter.fit(potential, {equilateralTriangle}, rules);
    const auto& coeffs = potential.getComponents()[0]->getCoefficients();
    ASSERT_NEAR(coeffs[0], 0.15381640115291151, 1e-8);
}

Atoms setupPythagorean() {
    return Atoms({{0.0, 0.0, 0.0}, {4.0, 0.0, 0.0}, {0.0, 3.0, 0.0}}, {Species("Fe"), Species("Fe"), Species("Fe")},
                 Lattice{{20.0, 0.0, 0.0}, {0.0, 20.0, 0.0}, {0.0, 0.0, 20.0}});
}

TEST(TestQRGapFit, pythagorian3bQuipCompatibility) {
    auto pythagorian = setupPythagorean();
    pythagorian.setEnergy(1.0);
    pythagorian.setForces({{1.0, 0.0, 0.0}, {0.0, -1.0, 0.0}, {0.5, 0.5, 0.0}});
    auto potential = create3bPotential(1.0, 1.0, 0.6, 10.0, {{6.0, 0.0, 3.0}});
    auto rules = SimpleRegularizationRules(1.0, 1.0, 1.0, 1.0);
    QRGapFit fitter;
    fitter.fit(potential, {pythagorian}, rules);
    const auto& coeffs = potential.getComponents()[0]->getCoefficients();
    ASSERT_NEAR(coeffs[0], -0.13044210897182607, 1e-8);
}
