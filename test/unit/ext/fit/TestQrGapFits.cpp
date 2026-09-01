#include <gtest/gtest.h>

#include "../../../../src/jgap/core/io/xyz/XYZData.hpp"
#include "jgap/core/atomic/Atoms.hpp"
#include "jgap/core/atomic/neighbours/NeighbourLists.hpp"
#include "jgap/core/atomic/species/composition/Species2Atomic.hpp"
#include "jgap/core/atomic/species/composition/Species2Sorted.hpp"
#include "jgap/core/atomic/species/composition/Species3AtomicSorted.hpp"
#include "jgap/core/cutoff/CosCutoff.hpp"
#include "jgap/core/fit/gap/regularization/SimpleRegularizationRules.hpp"
#include "jgap/core/kernels/SquaredExpKernel.hpp"
#include "jgap/core/potentials/gap/GapPotential.hpp"
#include "jgap/core/potentials/gap/component/ManyBodyGapComponent.hpp"
#include "jgap/core/potentials/gap/component/ThreeBodyGapComponent.hpp"
#include "jgap/core/potentials/gap/component/TwoBodyGapComponent.hpp"
#include "jgap/core/sparsification/HistogramUniformSparsifier.hpp"
#include "jgap/core/transform/manybody/TwoBodySum.hpp"
#include "jgap/core/transform/nbody/2b/PairDistanceTransformation.hpp"
#include "jgap/core/transform/nbody/2b/eam/PolycutoffPairFunction.hpp"
#include "jgap/core/transform/nbody/3b/Angle3bTransformation.hpp"
#include "jgap/experimental/fit/gap/BlockIncrementalQRGapFit.hpp"
#include "jgap/experimental/fit/gap/ElementalQRGapFit.hpp"
#include "jgap/ext/fit/gap/QRGapFit.hpp"
#include "jgap/utils/gap/GapComponentUtils.hpp"

using namespace jgap;

template<typename T>
class QrGapFits : public ::testing::Test {
protected:
    T fitter;
};

struct TestBlockIncrementalQRGapFit : public BlockIncrementalQRGapFit {
    TestBlockIncrementalQRGapFit() : BlockIncrementalQRGapFit(1e-8, 1.0) {}
};

using QrFitterTypes = ::testing::Types<QRGapFit, TestBlockIncrementalQRGapFit>;
TYPED_TEST_SUITE(QrGapFits, QrFitterTypes);

Atoms setupEquilateralTriangle() {
    return Atoms(
        {{0.0, 0.0, 0.0}, {3.0, 0.0, 0.0}, {1.5, 2.598, 0.0}},
        {Species("Fe"), Species("Fe"), Species("Fe")},
        Lattice{{100.0, 0.0, 0.0}, {0.0, 100.0, 0.0}, {0.0, 0.0, 100.0}}
    );
}

TYPED_TEST(QrGapFits, twoBodyEquilateralTriangleAtEquilibriumQuipCompatibility) {
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

    TypeParam fitter;
    std::vector<Atoms> train_data = {equilateralTriangle};
    fitter.fit(potential, train_data, regularization_rules.determineForAll(train_data));

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

    auto component = TwoBodyGapComponent<2, SquaredExpKernel<1, 1>>(
        Species2Sorted("Fe", "Fe"), std::move(trans), kernel, sparse_points
    );

    return GapPotential({component});
}

Atoms setupTwoAtoms() {
    return Atoms(
        {{0.0, 0.0, 0.0}, {3.0, 0.0, 0.0}},
        {Species("Fe"), Species("Fe")},
        Lattice{{100.0, 0.0, 0.0}, {0.0, 100.0, 0.0}, {0.0, 0.0, 100.0}}
    );
}

TYPED_TEST(QrGapFits, twoAtomsWithForceQuipCompatibility1) {
    auto twoAtoms = setupTwoAtoms();
    twoAtoms.setEnergy(1.0);
    twoAtoms.setForces({{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}});
    auto potential = create2bPotential({2.0, 4.0});
    auto rules = SimpleRegularizationRules(1.0, 1.0, 5.0, 5.0);
    TypeParam fitter;
    std::vector<Atoms> train_data = {twoAtoms};
    fitter.fit(potential, train_data, rules.determineForAll(train_data));
    const auto& coeffs = potential.getComponents()[0]->getCoefficients();
    ASSERT_NEAR(coeffs[0], 0.23266775926220731, 1e-6);
    ASSERT_NEAR(coeffs[1], 0.23266775926220731, 1e-6);
}

TYPED_TEST(QrGapFits, twoAtomsWithForceQuipCompatibility2) {
    auto twoAtoms = setupTwoAtoms();
    twoAtoms.setEnergy(0.0);
    twoAtoms.setForces({{1.0, 0.0, 0.0}, {0.0, 0.0, 0.0}});
    auto potential = create2bPotential({2.0});
    auto rules = SimpleRegularizationRules(10000000.0, 1.0, 5.0, 5.0);
    TypeParam fitter;
    std::vector<Atoms> train_data = {twoAtoms};
    fitter.fit(potential, train_data, rules.determineForAll(train_data));
    const auto& coeffs = potential.getComponents()[0]->getCoefficients();
    ASSERT_NEAR(coeffs[0], -0.30764655994411466, 1e-6);
}

TYPED_TEST(QrGapFits, twoAtomsWithForceQuipCompatibility3) {
    Atoms twoAtoms(
        {{0.0, 0.0, 0.0}, {1.5, 2.598, 0.0}},
        {Species("Fe"), Species("Fe")},
        Lattice{{100.0, 0.0, 0.0}, {0.0, 100.0, 0.0}, {0.0, 0.0, 100.0}}
    );
    twoAtoms.setEnergy(0.0);
    twoAtoms.setForces({{1.0, 0.0, 0.0}, {0.0, 0.0, 0.0}});
    auto potential = create2bPotential({2.0});
    auto rules = SimpleRegularizationRules(10000000.0, 1.0, 5.0, 5.0);
    TypeParam fitter;
    std::vector<Atoms> train_data = {twoAtoms};
    fitter.fit(potential, train_data, rules.determineForAll(train_data));
    const auto& coeffs = potential.getComponents()[0]->getCoefficients();
    ASSERT_NEAR(coeffs[0], -0.15382666452610533, 1e-6);
}

TYPED_TEST(QrGapFits, twoAtomsWithForceQuipCompatibility4) {
    Atoms twoAtoms(
        {{0.0, 0.0, 0.0}, {1.5, 2.598, 0.0}},
        {Species("Fe"), Species("Fe")},
        Lattice{{100.0, 0.0, 0.0}, {0.0, 100.0, 0.0}, {0.0, 0.0, 100.0}}
    );
    twoAtoms.setEnergy(0.0);
    twoAtoms.setForces({{1.0, 1.0, 1.0}, {-1.0, -1.0, -1.0}});
    auto potential = create2bPotential({2.0});
    auto rules = SimpleRegularizationRules(10000000.0, 2.0, 5.0, 5.0);
    TypeParam fitter;
    std::vector<Atoms> train_data = {twoAtoms};
    fitter.fit(potential, train_data, rules.determineForAll(train_data));
    const auto& coeffs = potential.getComponents()[0]->getCoefficients();
    ASSERT_NEAR(coeffs[0], -0.47733536744854282, 1e-6);
}

TYPED_TEST(QrGapFits, twoAtomsWithForceQuipCompatibility5) {
    Atoms twoAtoms(
        {{0.0, 0.0, 0.0}, {1.5, 2.598, 0.0}},
        {Species("Fe"), Species("Fe")},
        Lattice{{100.0, 0.0, 0.0}, {0.0, 100.0, 0.0}, {0.0, 0.0, 100.0}}
    );
    twoAtoms.setEnergy(0.0);
    twoAtoms.setForces({{1.0, 1.0, 1.0}, {-1.0, -1.0, -1.0}});
    auto potential = create2bPotential({2.0, 2.5, 4.0});
    auto rules = SimpleRegularizationRules(1.0, 2.0, 5.0, 5.0);
    TypeParam fitter;
    std::vector<Atoms> train_data = {twoAtoms};
    fitter.fit(potential, train_data, rules.determineForAll(train_data));
    const auto& coeffs = potential.getComponents()[0]->getCoefficients();
    ASSERT_NEAR(coeffs[0], -0.24567921601436607, 1e-6);
    ASSERT_NEAR(coeffs[1], -0.10998748761736651, 1e-6);
    ASSERT_NEAR(coeffs[2], 0.38697171383300527, 1e-6);
}

GapPotential createEamPotential(Real theta, Real delta, Real r_min, Real cutoff, const std::vector<Real>& sparse_pts) {
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

TYPED_TEST(QrGapFits, twoAtomsEamQuipCompatibility) {
    auto twoAtoms = setupTwoAtoms();
    twoAtoms.setEnergy(1.0);
    twoAtoms.setForces({{1.0, 0.0, 0.0}, {-1.0, 0.0, 0.0}});
    auto potential = createEamPotential(1.0, 1.0, 0.0, 5.0, {1.0});
    auto rules = SimpleRegularizationRules(1.0, 1.0, 1.0, 1.0);
    TypeParam fitter;
    std::vector<Atoms> train_data = {twoAtoms};
    fitter.fit(potential, train_data, rules.determineForAll(train_data));
    const auto& coeffs = potential.getComponents()[0]->getCoefficients();
    ASSERT_NEAR(coeffs[0], 0.017637586136984833, 1e-6);
}

TYPED_TEST(QrGapFits, eamQuipCompatibilityRealBox) {
    auto box = Atoms::readAtoms("test/resources/xyz-samples/fe-only.xyz")[15];
    box.eraseVirials();
    auto potential = createEamPotential(3.0, 2.0, 0.0, 5.0, {1.0, 2.5, 4.0});
    auto rules = SimpleRegularizationRules(0.001, 0.05, 1.0, 1.0);
    TypeParam fitter;
    std::vector<Atoms> train_data = {box};
    fitter.fit(potential, train_data, rules.determineForAll(train_data));
    const auto& coeffs = potential.getComponents()[0]->getCoefficients();

    ASSERT_NEAR(coeffs[0], -5.0619173053329485, 1e-6);
    ASSERT_NEAR(coeffs[1], -4.4169841885357615, 1e-6);
    ASSERT_NEAR(coeffs[2], -0.16117464517539212, 1e-6);
}

GapPotential create3bPotential(
    Real theta, Real delta, Real cutoff_transition_width, Real cutoff, const std::vector<Vector3>& sparsePts
) {
    auto trans = Angle3bTransformation(CosCutoff(cutoff, cutoff_transition_width));
    auto kernel = SquaredExpKernel<3, 1>(delta, std::array<Real, 3>{theta, theta, theta});
    std::vector<Descriptor<4>> sparse_points;
    for (const auto& q: sparsePts) {
        sparse_points.push_back({q.x, q.y, q.z, 1.0});
    }
    auto component = ThreeBodyGapComponent<4, SquaredExpKernel<3, 1>>(
        Species3AtomicSorted("Fe", "Fe", "Fe"), std::move(trans), std::move(kernel), sparse_points
    );

    return GapPotential({component});
}

TYPED_TEST(QrGapFits, equilateralTriangle3bQuipCompatibility) {
    auto equilateralTriangle = setupEquilateralTriangle();
    equilateralTriangle.setEnergy(1.0);
    equilateralTriangle.setForces({{1.0, 0.0, 0.0}, {0.0, -1.0, 0.0}, {0.5, 0.5, 0.0}});
    auto potential = create3bPotential(1.0, 1.0, 0.6, 10.0, {{6.0, 0.0, 3.0}});
    auto rules = SimpleRegularizationRules(1.0, 1.0, 1.0, 1.0);
    TypeParam fitter;
    std::vector<Atoms> train_data = {equilateralTriangle};
    fitter.fit(potential, train_data, rules.determineForAll(train_data));
    const auto& coeffs = potential.getComponents()[0]->getCoefficients();
    ASSERT_NEAR(coeffs[0], 0.15381640115291151, 1e-8);
}

Atoms setupPythagorean() {
    return Atoms(
        {{0.0, 0.0, 0.0}, {4.0, 0.0, 0.0}, {0.0, 3.0, 0.0}},
        {Species("Fe"), Species("Fe"), Species("Fe")},
        Lattice{{20.0, 0.0, 0.0}, {0.0, 20.0, 0.0}, {0.0, 0.0, 20.0}}
    );
}

TYPED_TEST(QrGapFits, pythagorian3bQuipCompatibility) {
    auto pythagorian = setupPythagorean();
    pythagorian.setEnergy(1.0);
    pythagorian.setForces({{1.0, 0.0, 0.0}, {0.0, -1.0, 0.0}, {0.5, 0.5, 0.0}});
    auto potential = create3bPotential(1.0, 1.0, 0.6, 10.0, {{6.0, 0.0, 3.0}});
    auto rules = SimpleRegularizationRules(1.0, 1.0, 1.0, 1.0);
    TypeParam fitter;
    std::vector<Atoms> train_data = {pythagorian};
    fitter.fit(potential, train_data, rules.determineForAll(train_data));
    const auto& coeffs = potential.getComponents()[0]->getCoefficients();
    ASSERT_NEAR(coeffs[0], -0.13044210897182607, 1e-8);
}

TYPED_TEST(QrGapFits, MissingRegularizationThrows) {
    auto twoAtoms = setupTwoAtoms();
    twoAtoms.setEnergy(1.0);
    auto potential = create2bPotential({2.0});
    TypeParam fitter;

    // Empty sigmas
    EXPECT_THROW(fitter.fit(potential, {twoAtoms}, {}), std::runtime_error);

    // Missing energy sigma
    Regularization empty_reg;
    EXPECT_THROW(fitter.fit(potential, {twoAtoms}, {empty_reg}), std::runtime_error);
}

TEST(SplitQRGapFitTest, NonZeroCovarianceForSpecies) {
    auto trans2 = PairDistanceTransformation(CosCutoff(5.0, 1.0));
    auto kernel2 = SquaredExpKernel<1, 1>(1.0, std::array<Real, 1>{1.0});
    TwoBodyGapComponent<2, SquaredExpKernel<1, 1>> comp2(
        Species2Sorted("Fe", "Ni"), trans2, kernel2, std::vector<Descriptor<2>>{{2.5, 1.0}}
    );
    EXPECT_EQ(comp2.nonZeroCovarianceFor(), (std::set<Species>{"Fe", "Ni"}));

    auto trans3 = Angle3bTransformation(CosCutoff(5.0, 1.0));
    auto kernel3 = SquaredExpKernel<3, 1>(1.0, std::array<Real, 3>{1.0, 1.0, 1.0});
    ThreeBodyGapComponent<4, SquaredExpKernel<3, 1>> comp3(
        Species3AtomicSorted("Fe", "Ni", "Cu"), trans3, kernel3, std::vector<Descriptor<4>>{{2.5, 2.5, 2.5, 1.0}}
    );
    EXPECT_EQ(comp3.nonZeroCovarianceFor(), (std::set<Species>{"Fe", "Ni", "Cu"}));

    TwoBodySum<1> sum(Species("Fe"));
    sum.extend(Species2Atomic("Fe", "Ni"), PolycutoffPairFunction(5.0, 0.0));
    auto kernel_eam = SquaredExpKernel<1, 0>(1.0, std::array<Real, 1>{1.0});
    ManyBodyGapComponent<1, SquaredExpKernel<1, 0>> comp_mb(
        sum, kernel_eam, std::vector<Descriptor<1>>{{1.0}}
    );
    EXPECT_EQ(comp_mb.nonZeroCovarianceFor(), (std::set<Species>{"Fe"}));
}

TEST(ElementalQRGapFitTest, MultiSpeciesSplitFitMatchesQRGapFit) {
    // 1. Fe dimer
    Atoms fe_dimer(
        {{0.0, 0.0, 0.0}, {2.5, 0.0, 0.0}},
        {Species("Fe"), Species("Fe")},
        Lattice{{20.0, 0.0, 0.0}, {0.0, 20.0, 0.0}, {0.0, 0.0, 20.0}}
    );
    fe_dimer.setEnergy(-1.5);
    fe_dimer.setForces({{0.1, 0.0, 0.0}, {-0.1, 0.0, 0.0}});

    // 2. Ni dimer
    Atoms ni_dimer(
        {{0.0, 0.0, 0.0}, {2.4, 0.0, 0.0}},
        {Species("Ni"), Species("Ni")},
        Lattice{{20.0, 0.0, 0.0}, {0.0, 20.0, 0.0}, {0.0, 0.0, 20.0}}
    );
    ni_dimer.setEnergy(-2.0);
    ni_dimer.setForces({{0.2, 0.0, 0.0}, {-0.2, 0.0, 0.0}});

    // 3. Fe-Ni mixed dimer
    Atoms feni_dimer(
        {{0.0, 0.0, 0.0}, {2.45, 0.0, 0.0}},
        {Species("Fe"), Species("Ni")},
        Lattice{{20.0, 0.0, 0.0}, {0.0, 20.0, 0.0}, {0.0, 0.0, 20.0}}
    );
    feni_dimer.setEnergy(-1.8);
    feni_dimer.setForces({{0.05, 0.0, 0.0}, {-0.05, 0.0, 0.0}});

    std::vector<Atoms> train_data = {fe_dimer, ni_dimer, feni_dimer};

    auto trans = PairDistanceTransformation(CosCutoff(6.0, 1.0));
    auto kernel = SquaredExpKernel<1, 1>(1.0, std::array<Real, 1>{1.0});

    TwoBodyGapComponent<2, SquaredExpKernel<1, 1>> comp_fe(
        Species2Sorted("Fe", "Fe"), trans, kernel, std::vector<Descriptor<2>>{{2.5, 1.0}}
    );
    TwoBodyGapComponent<2, SquaredExpKernel<1, 1>> comp_ni(
        Species2Sorted("Ni", "Ni"), trans, kernel, std::vector<Descriptor<2>>{{2.4, 1.0}}
    );
    TwoBodyGapComponent<2, SquaredExpKernel<1, 1>> comp_feni(
        Species2Sorted("Fe", "Ni"), trans, kernel, std::vector<Descriptor<2>>{{2.45, 1.0}}
    );

    GapPotential pot_qr({comp_fe, comp_ni, comp_feni});
    GapPotential pot_split({comp_fe, comp_ni, comp_feni});

    SimpleRegularizationRules rules(1.0, 1.0, 1.0, 1.0);
    auto sigmas = rules.determineForAll(train_data);

    QRGapFit qr_fitter(1e-8);
    qr_fitter.fit(pot_qr, train_data, sigmas);

    ElementalQRGapFit split_fitter(1e-8, 1.0);
    split_fitter.fit(pot_split, train_data, sigmas);

    // Verify coefficients match
    ASSERT_EQ(pot_qr.getComponents().size(), pot_split.getComponents().size());
    for (size_t c = 0; c < pot_qr.getComponents().size(); ++c) {
        const auto& c_qr = pot_qr.getComponents()[c]->getCoefficients();
        const auto& c_sp = pot_split.getComponents()[c]->getCoefficients();
        ASSERT_EQ(c_qr.size(), c_sp.size());
        for (size_t i = 0; i < c_qr.size(); ++i) {
            EXPECT_NEAR(c_qr[i], c_sp[i], 1e-10);
        }
    }

    // Verify both potentials evaluate energies and forces identically on all frames
    for (const auto& atoms: train_data) {
        auto e_qr = pot_qr.calculateEnergy(atoms);
        auto e_split = pot_split.calculateEnergy(atoms);
        EXPECT_NEAR(e_qr.value, e_split.value, 1e-10);
        for (size_t a = 0; a < atoms.nAtoms(); ++a) {
            EXPECT_NEAR(e_qr.forces[a].x, e_split.forces[a].x, 1e-10);
            EXPECT_NEAR(e_qr.forces[a].y, e_split.forces[a].y, 1e-10);
            EXPECT_NEAR(e_qr.forces[a].z, e_split.forces[a].z, 1e-10);
        }
    }
}

TEST(ElementalQRGapFitTest, FeNiTrainDatasetConsistencyAcrossFitters) {
    auto all_atoms = Atoms::readAtoms("test/resources/xyz-samples/feni-train.xyz");
    std::vector<Atoms> train_data;
    std::vector<Atoms> fe_only, ni_only, feni_both;

    for (const auto& a: all_atoms) {
        std::set<Species> sp(a.getSpecies().begin(), a.getSpecies().end());
        if (sp == std::set<Species>{Species("Fe")} && fe_only.size() < 10) {
            fe_only.push_back(a);
        } else if (sp == std::set<Species>{Species("Ni")} && ni_only.size() < 10) {
            ni_only.push_back(a);
        } else if (sp == std::set<Species>{Species("Fe"), Species("Ni")} && feni_both.size() < 10) {
            feni_both.push_back(a);
        }
    }
    ASSERT_EQ(fe_only.size(), 10);
    ASSERT_EQ(ni_only.size(), 10);
    ASSERT_EQ(feni_both.size(), 10);

    train_data.insert(train_data.end(), fe_only.begin(), fe_only.end());
    train_data.insert(train_data.end(), ni_only.begin(), ni_only.end());
    train_data.insert(train_data.end(), feni_both.begin(), feni_both.end());

    // 4 sparse pts for 2b
    auto trans2 = PairDistanceTransformation(CosCutoff(4.5, 1.0));
    auto kernel2 = SquaredExpKernel<1, 1>(10.0, {1.0});
    auto sparsifier2 = HistogramUniformSparsifier<2>(42, 4, std::array{true, false});
    auto comps2 = utils::createTwoBodyComponents<2, SquaredExpKernel<1, 1>>(train_data, trans2, kernel2, sparsifier2);

    // 4 sparse pts for eam
    auto eam_pf = PolycutoffPairFunction(4.5, 0.0);
    auto kernel_eam = SquaredExpKernel<1, 0>(1.0, {1.0});
    auto sparsifier_eam = HistogramUniformSparsifier<1>(42, 4, std::nullopt, std::nullopt, Descriptor<1>{0.05});
    auto comps_eam = utils::createEamComponents<SquaredExpKernel<1, 0>>(
        eam_pf, kernel_eam, sparsifier_eam, train_data, EamMode::Blind
    );

    // 10 sparse pts for 3b
    auto trans3 = Angle3bTransformation(CosCutoff(4.0, 1.0));
    auto kernel3 = SquaredExpKernel<3, 1>(1.0, {1.0, 1.0, 1.0});
    auto sparsifier3 = HistogramUniformSparsifier<4>(42, 10, std::array{true, true, true, false});
    auto comps3 = utils::createThreeBodyComponents<4, SquaredExpKernel<3, 1>>(train_data, trans3, kernel3, sparsifier3);

    std::vector<ValuePtr<GapComponent>> all_comps;
    for (const auto& c: comps2) all_comps.emplace_back(c);
    for (const auto& c: comps_eam) all_comps.emplace_back(c);
    for (const auto& c: comps3) all_comps.emplace_back(c);

    GapPotential pot_qr;
    pot_qr.components = all_comps;

    GapPotential pot_stream;
    for (const auto& c: all_comps) pot_stream.components.push_back(c);

    GapPotential pot_split;
    for (const auto& c: all_comps) pot_split.components.push_back(c);

    SimpleRegularizationRules rules(1.0, 1.0, 1.0, 1.0);
    auto sigmas = rules.determineForAll(train_data);

    // Use small RAM limit so that chunks/splits are exercised (e.g. 0.0002 GB ~ 200 KB)
    const double test_ram_gb = 0.0002;

    QRGapFit qr_fitter(1e-8);
    qr_fitter.fit(pot_qr, train_data, sigmas);

    BlockIncrementalQRGapFit stream_fitter(1e-8, test_ram_gb);
    stream_fitter.fit(pot_stream, train_data, sigmas);

    ElementalQRGapFit split_fitter(1e-8, test_ram_gb);
    split_fitter.fit(pot_split, train_data, sigmas);

    // 1. Verify pot_qr and pot_stream match component-wise directly
    for (size_t c = 0; c < pot_qr.getComponents().size(); ++c) {
        const auto& c_qr = pot_qr.getComponents()[c]->getCoefficients();
        const auto& c_st = pot_stream.getComponents()[c]->getCoefficients();
        ASSERT_EQ(c_qr.size(), c_st.size());
        for (size_t i = 0; i < c_qr.size(); ++i) {
            EXPECT_NEAR(c_qr[i], c_st[i], 1e-8);
        }
    }

    // 2. Verify all three potentials evaluate to identical energies and forces on all 30 structures
    for (const auto& atoms: train_data) {
        auto e_qr = pot_qr.calculateEnergy(atoms);
        auto e_st = pot_stream.calculateEnergy(atoms);
        auto e_sp = pot_split.calculateEnergy(atoms);

        EXPECT_NEAR(e_qr.value, e_st.value, 1e-7);
        EXPECT_NEAR(e_qr.value, e_sp.value, 1e-7);

        for (size_t a = 0; a < atoms.nAtoms(); ++a) {
            EXPECT_NEAR(e_qr.forces[a].x, e_st.forces[a].x, 1e-7);
            EXPECT_NEAR(e_qr.forces[a].y, e_st.forces[a].y, 1e-7);
            EXPECT_NEAR(e_qr.forces[a].z, e_st.forces[a].z, 1e-7);

            EXPECT_NEAR(e_qr.forces[a].x, e_sp.forces[a].x, 1e-7);
            EXPECT_NEAR(e_qr.forces[a].y, e_sp.forces[a].y, 1e-7);
            EXPECT_NEAR(e_qr.forces[a].z, e_sp.forces[a].z, 1e-7);
        }
    }
}


