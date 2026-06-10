#include <gtest/gtest.h>

#include "core/fit/gap/QRGapFit.hpp"
#include "core/atomic/Atoms.hpp"
#include "core/atomic/io/XYZData.hpp"
#include "core/atomic/neighbours/NeighbourList.hpp"
#include "core/potentials/gap/component/NBodyGapComponent.hpp"
#include "core/potentials/gap/component/ManyBodyGapComponent.hpp"
#include "core/potentials/gap/GapPotential.hpp"
#include "core/transform/2b/TwoBodyTransformation.hpp"
#include "core/transform/3b/Angle3bTransformation.hpp"
#include "core/transform/eam/PolycutoffPairFunction.hpp"
#include "../../../src/core/transform/aggregated/TransformationAggregator.hpp"
#include "core/kernels/SquaredExpKernel.hpp"
#include "core/cutoff/CosCutoff.hpp"
#include "core/fit/gap/regularization/SimpleRegularizationRules.hpp"
#include "core/transform/aggregated/TransformationAggregatorImpl.hpp"

using namespace jgap;

Atoms setupEquilateralTriangle() {
    return Atoms(
        { {0.0, 0.0, 0.0}, {3.0, 0.0, 0.0}, {1.5, 2.598, 0.0} },
        { Species("Fe"), Species("Fe"), Species("Fe") },
        Lattice{ {100.0, 0.0, 0.0}, {0.0, 100.0, 0.0}, {0.0, 0.0, 100.0} }
    );
}

TEST(TestQRGapFit, twoBodyEquilateralTriangleAtEquilibriumQuipCompatibility) {
    auto equilateralTriangle = setupEquilateralTriangle();
    equilateralTriangle.setEnergy(1.0);
    equilateralTriangle.setForces({{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}});

    auto trans = std::make_unique<TwoBodyTransformation>(std::make_unique<CosCutoff>(10.0, 0.7));
    auto kernel = std::make_unique<SquaredExpKernel<1, 1>>(1.0, std::array<Real, 1>{1.0});
    std::vector<Descriptor<2>> sparse_points = { {{{3.0, 1.0}}} };
    auto component = std::make_unique<NBodyGapComponent<2, 2, Symmetric>>(
        SpeciesSet<2, Symmetric>{"Fe", "Fe"}, std::move(trans), std::move(kernel), sparse_points
        );

    std::vector<GapComponent::Ptr> components;
    components.push_back(std::move(component));
    auto potential = GapPotential(std::move(components));
    auto regularization_rules = SimpleRegularizationRules(3.0, 10.0, 5.0, 5.0);

    QRGapFit fitter;
    fitter.fit(potential, {equilateralTriangle}, regularization_rules);

    const auto& fitted_coeffs = potential.getComponents()[0]->getCoefficients();
    double c = fitted_coeffs[0];
    auto totalVariance = 6.0 / c - 36.0;
    ASSERT_NEAR(totalVariance, 27, 1e-6);
}

GapPotential create2bPotential(const std::vector<Real>& sparsePts) {
    auto trans = std::make_unique<TwoBodyTransformation>(std::make_unique<CosCutoff>(10.0, 0.7));
    auto kernel = std::make_unique<SquaredExpKernel<1, 1>>(1.0, std::array<Real, 1>{1.0});

    std::vector<Descriptor<2>> sparse_points;
    CosCutoff ref_cutoff(10.0, 0.7);
    for (Real r : sparsePts) {
        sparse_points.push_back({{{r, ref_cutoff.evaluate(r)}}});
    }

    auto component = std::make_unique<NBodyGapComponent<2, 2, Symmetric>>(
        SpeciesSet<2, Symmetric>{"Fe", "Fe"}, std::move(trans), std::move(kernel), sparse_points);
    std::vector<GapComponent::Ptr> components;
    components.push_back(std::move(component));
    return GapPotential{std::move(components)};
}

Atoms setupTwoAtoms() {
    return Atoms(
        { {0.0, 0.0, 0.0}, {3.0, 0.0, 0.0} },
        { Species("Fe"), Species("Fe") },
        Lattice{ {100.0, 0.0, 0.0}, {0.0, 100.0, 0.0}, {0.0, 0.0, 100.0} }
    );
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
    Atoms twoAtoms({ {0.0, 0.0, 0.0}, {1.5, 2.598, 0.0} }, { Species("Fe"), Species("Fe") }, Lattice{ {100,0,0}, {0,100,0}, {0,0,100} });
    twoAtoms.setEnergy(0.0);
    twoAtoms.setForces({{1.0, 0.0, 0.0}, {0.0, 0.0, 0.0}});
    auto potential = create2bPotential({2.0});
    auto rules = SimpleRegularizationRules(10000000.0, 1.0, 5.0, 5.0);
    QRGapFit fitter;
    fitter.fit(potential, {twoAtoms}, rules);
    const auto& coeffs = potential.getComponents()[0]->getCoefficients();
    ASSERT_NEAR(coeffs[0], -.15382666452610533, 1e-6);
}

TEST(TestQRGapFit, twoAtomsWithForceQuipCompatibility4) {
    Atoms twoAtoms({ {0.0, 0.0, 0.0}, {1.5, 2.598, 0.0} }, { Species("Fe"), Species("Fe") }, Lattice{ {100,0,0}, {0,100,0}, {0,0,100} });
    twoAtoms.setEnergy(0.0);
    twoAtoms.setForces({{1.0, 1.0, 1.0}, {-1.0, -1.0, -1.0}});
    auto potential = create2bPotential({2.0});
    auto rules = SimpleRegularizationRules(10000000.0, 2.0, 5.0, 5.0);
    QRGapFit fitter;
    fitter.fit(potential, {twoAtoms}, rules);
    const auto& coeffs = potential.getComponents()[0]->getCoefficients();
    ASSERT_NEAR(coeffs[0], -.47733536744854282, 1e-6);
}

TEST(TestQRGapFit, twoAtomsWithForceQuipCompatibility5) {
    Atoms twoAtoms({ {0.0, 0.0, 0.0}, {1.5, 2.598, 0.0} }, { Species("Fe"), Species("Fe") }, Lattice{ {100,0,0}, {0,100,0}, {0,0,100} });
    twoAtoms.setEnergy(0.0);
    twoAtoms.setForces({{1.0, 1.0, 1.0}, {-1.0, -1.0, -1.0}});
    auto potential = create2bPotential({2.0, 2.5, 4.0});
    auto rules = SimpleRegularizationRules(1.0, 2.0, 5.0, 5.0);
    QRGapFit fitter;
    fitter.fit(potential, {twoAtoms}, rules);
    const auto& coeffs = potential.getComponents()[0]->getCoefficients();
    ASSERT_NEAR(coeffs[0], -.24567921601436607, 1e-6);
    ASSERT_NEAR(coeffs[1], -.10998748761736651, 1e-6);
    ASSERT_NEAR(coeffs[2], .38697171383300527, 1e-6);
}

GapPotential createEamPotential(double theta, double delta, double r_min, double cutoff, const std::vector<Real>& sparse_pts) {
    auto trans = PolycutoffPairFunction(cutoff, r_min, 1.0);

    auto eam_aggregator = std::make_unique<TransformationAggregatorImpl<1, 2>>(Species("Fe"));

    eam_aggregator->extend(SpeciesSet<2, HasCentralAtom>{"Fe", "Fe"}, std::make_unique<PolycutoffPairFunction>(trans));
    auto kernel = std::make_unique<SquaredExpKernel<1, 0>>(delta, std::array<Real, 1>{theta});
    std::vector<Descriptor<1>> sparse_points;
    for (Real r : sparse_pts) {
        sparse_points.push_back({{{r}}});
    }
    auto component = std::make_unique<ManyBodyGapComponent<1>>(
        std::move(eam_aggregator), std::move(kernel), sparse_points
    );
    std::vector<GapComponent::Ptr> components;
    components.push_back(std::move(component));
    return GapPotential(std::move(components));
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
    ASSERT_NEAR(coeffs[0], .17637586136984833E-001, 1e-6);
}

TEST(TestQRGapFit, eamQuipCompatibilityRealBox) {
    auto box = readAtoms("test/resources/xyz-samples/fe-only.xyz")[15];
    box.properties.erase("virial"); // initial tests ignored virials unfortunately
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
    ASSERT_NEAR(coeffs[2], -.16117464517539212, 1e-6);
}

GapPotential create3bPotential(double theta, double delta, double cutoff_transition_width, double cutoff, const std::vector<Vector3>& sparsePts) {
    auto trans = std::make_unique<Angle3bTransformation>(std::make_unique<CosCutoff>(cutoff, cutoff_transition_width));
    auto kernel = std::make_unique<SquaredExpKernel<3, 1>>(delta, std::array<Real, 3>{theta, theta, theta});
    std::vector<Descriptor<4>> sparse_points;
    for (const auto& q : sparsePts) {
        sparse_points.push_back({{{q.x, q.y, q.z, 1.0}}}); // Assume f_cut_prod=1 for sparse points
    }
    auto component = std::make_unique<NBodyGapComponent<4, 3, HasCentralAtom>>(
        SpeciesSet<3, HasCentralAtom>{"Fe", "Fe", "Fe"},
        std::move(trans),
        std::move(kernel),
        sparse_points
        );
    std::vector<GapComponent::Ptr> components;
    components.push_back(std::move(component));
    return GapPotential(std::move(components));
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
    ASSERT_NEAR(coeffs[0], .15381640115291151, 1e-8);
}

Atoms setupPythagorean() {
    return Atoms(
        { {0,0,0}, {4,0,0}, {0,3,0} },
        { Species("Fe"), Species("Fe"), Species("Fe") },
        Lattice{ {20,0,0}, {0,20,0}, {0,0,20} }
    );
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
    ASSERT_NEAR(coeffs[0], -.13044210897182607, 1e-8);
}

