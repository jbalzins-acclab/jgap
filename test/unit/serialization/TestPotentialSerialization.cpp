#include <gtest/gtest.h>

#include "jgap/core/atomic/Atoms.hpp"
#include "jgap/core/cutoff/CosCutoff.hpp"
#include "jgap/core/kernels/SquaredExpKernel.hpp"
#include "jgap/core/potentials/CompositePotential.hpp"
#include "jgap/core/potentials/gap/GapPotential.hpp"
#include "jgap/core/potentials/gap/component/TwoBodyGapComponent.hpp"
#include "jgap/core/potentials/isolated/IsolatedAtomPotential.hpp"
#include "jgap/core/potentials/spline/SplinePairPotential.hpp"
#include "jgap/core/potentials/tabgap/TabGapPotential.hpp"
#include "jgap/core/potentials/zbl/ZblPotential.hpp"
#include "jgap/core/transform/nbody/2b/PairDistanceTransformation.hpp"
#include "jgap/serialization/SerializationRegistry.hpp"

using namespace jgap;

namespace {
    std::string tmpFile(const std::string& name) { return testing::TempDir() + "/" + name; }

    template<typename Obj>
    ValuePtr<Potential> roundTrip(const Obj& potential, const std::string& file) {
        SerializationRegistry<Potential>::serialize(ValuePtr<Potential>(potential), file);
        return SerializationRegistry<Potential>::deserialize(file);
    }

    // Most potentials are compared by predicting on the same small structure before and after the
    // round-trip: identical energy/forces is strong evidence the whole object tree was restored.
    void expectSameEnergy(const Potential& a, const Potential& b, const Atoms& atoms) {
        const auto qa = a.calculateEnergy(atoms);
        const auto qb = b.calculateEnergy(atoms);
        EXPECT_NEAR(qa.value, qb.value, 1e-9);
        ASSERT_EQ(qa.forces.size(), qb.forces.size());
        for (size_t i = 0; i < qa.forces.size(); ++i) {
            EXPECT_NEAR(qa.forces[i].x, qb.forces[i].x, 1e-9);
            EXPECT_NEAR(qa.forces[i].y, qb.forces[i].y, 1e-9);
            EXPECT_NEAR(qa.forces[i].z, qb.forces[i].z, 1e-9);
        }
    }

    Atoms feNiPair() { return Atoms({{0, 0, 0}, {2.5, 0, 0}}, {Species("Fe"), Species("Ni")}); }
}

TEST(TestPotentialSerialization, IsolatedAtomPotential) {
    IsolatedAtomPotential original(std::map<Species, Real>{{Species("Fe"), -3.5}, {Species("Ni"), -2.1}});

    auto rt = roundTrip(original, tmpFile("isolated.h5"));
    auto restored = rt.as<IsolatedAtomPotential>();
    ASSERT_NE(restored, nullptr);

    const auto& energies = restored->getIsolatedEnergies();
    ASSERT_EQ(energies.size(), 2u);
    EXPECT_DOUBLE_EQ(energies.at(Species("Fe")), -3.5);
    EXPECT_DOUBLE_EQ(energies.at(Species("Ni")), -2.1);
}

TEST(TestPotentialSerialization, ZblPotential) {
    std::map<Species2Sorted, std::array<Real, 6>> coeffs = {{Species2Sorted("Fe,Ni"), {0.1, 1.2, 0.3, 1.4, 0.5, 1.6}}};
    ZblPotential original(coeffs, 8.0, 1.0);

    auto rt = roundTrip(original, tmpFile("zbl.h5"));
    auto restored = rt.as<ZblPotential>();
    ASSERT_NE(restored, nullptr);
    EXPECT_DOUBLE_EQ(restored->getCutoff(), 8.0);
    EXPECT_DOUBLE_EQ(restored->getCutoffTransitionWidth(), 1.0);

    const auto restored_coeffs = restored->getCoefficients();
    ASSERT_EQ(restored_coeffs.size(), 1u);
    const auto& c = restored_coeffs.at(Species2Sorted("Fe,Ni"));
    for (size_t i = 0; i < 6; ++i) {
        EXPECT_DOUBLE_EQ(c[i], coeffs.at(Species2Sorted("Fe,Ni"))[i]);
    }

    expectSameEnergy(original, *rt, feNiPair());
}

TEST(TestPotentialSerialization, SplinePairPotential) {
    SplinePairPotential original;
    std::vector<Real> r_vec = {1.0, 2.0, 3.0, 4.0};
    std::vector<Real> energies = {1.0, 0.25, 0.1, 0.0};
    original.extend(Species("Fe"), Species("Ni"), r_vec, energies);

    auto rt = roundTrip(original, tmpFile("spline_pair.h5"));
    auto restored = rt.as<SplinePairPotential>();
    ASSERT_NE(restored, nullptr);

    ASSERT_EQ(restored->getInterpolators().size(), 1u);
    const auto& spline = restored->getInterpolators().begin()->second;
    EXPECT_EQ(spline.getRVec(), r_vec);
    EXPECT_EQ(spline.getEnergies(), energies);

    expectSameEnergy(original, *rt, feNiPair());
}

TEST(TestPotentialSerialization, CompositePotential) {
    CompositePotential original{std::map<std::string, ValuePtr<Potential>>{
        {"isolated", IsolatedAtomPotential(std::map<Species, Real>{{Species("Fe"), -3.5}, {Species("Ni"), -2.1}})},
        {"zbl", ZblPotential(std::map<Species2Sorted, std::array<Real, 6>>{{Species2Sorted("Fe,Ni"),
                                                                            {0.1, 1.2, 0.3, 1.4, 0.5, 1.6}}},
                             8.0, 1.0)},
    }};

    auto rt = roundTrip(original, tmpFile("composite.h5"));
    auto restored = rt.as<CompositePotential>();
    ASSERT_NE(restored, nullptr);
    ASSERT_EQ(restored->getPotentials().size(), 2u);
    EXPECT_TRUE(restored->getPotentials().contains("isolated"));
    EXPECT_TRUE(restored->getPotentials().contains("zbl"));

    expectSameEnergy(original, *rt, feNiPair());
}

namespace {
    // A minimal fitted-shape GAP: one 2-body Fe-Fe component with explicit sparse points and
    // coefficients (no fitting needed) plus an isolated-atom external potential.
    GapPotential makeSmallGap() {
        std::vector<Descriptor<2>> sparse_points = {Descriptor<2>{{2.0, 1.0}}, Descriptor<2>{{2.6, 0.5}}};
        std::vector<Real> coefficients = {0.5, -0.3};

        GapPotential potential;
        auto transform_ptr = ValuePtr<TwoBodyTransformation<2>>(PairDistanceTransformation(CosCutoff(5.0, 1.0)));
        potential.addComponent(ValuePtr<GapComponent>(TwoBodyGapComponent<2, SquaredExpKernel<1, 1>>(
            Species2Sorted("Fe,Fe"), transform_ptr, SquaredExpKernel<1, 1>(10.0, {1.3}), sparse_points, coefficients)));
        potential.optional_external_potential = IsolatedAtomPotential(std::map<Species, Real>{{Species("Fe"), -3.5}});
        return potential;
    }

    Atoms fePair() { return Atoms({{0, 0, 0}, {2.5, 0, 0}}, {Species("Fe"), Species("Fe")}); }
}

TEST(TestPotentialSerialization, GapPotential) {
    GapPotential original = makeSmallGap();

    auto rt = roundTrip(original, tmpFile("gap.h5"));
    auto restored = rt.as<GapPotential>();
    ASSERT_NE(restored, nullptr);
    EXPECT_EQ(restored->getComponents().size(), original.getComponents().size());
    EXPECT_NE(restored->optional_external_potential.get(), nullptr);

    expectSameEnergy(original, *rt, fePair());
}

TEST(TestPotentialSerialization, TabGapPotential) {
    GapPotential gap = makeSmallGap();

    // Tabulate (the tabGAP path) then round-trip the tabulated potential through the special
    // jgap-flavoured tabGAP serializer (embedded .eam.fs, jgap=true marker).
    TabGapPotential original{gap.tabulate({
        .max_cutoffs = gap.getCutoffs(),
        .n_grid_2b = 200,
    })};

    auto rt = roundTrip(original, tmpFile("tabgap.h5"));
    auto restored = rt.as<TabGapPotential>();
    ASSERT_NE(restored, nullptr);

    expectSameEnergy(original, *rt, fePair());
}
