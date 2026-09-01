#include <gtest/gtest.h>
#include "jgap/core/atomic/Atoms.hpp"
#include "jgap/core/fit/gap/regularization/PerConfigTypeRegularizationRules.hpp"
#include "jgap/core/fit/gap/regularization/PerConfigTypeSigmas.hpp"
#include "jgap/core/fit/gap/regularization/ScaledRegularizationRules.hpp"
#include "jgap/core/fit/gap/regularization/SimpleRegularizationRules.hpp"

using namespace jgap;

TEST(TestScaledRegularizationRules, NoForcesDefaultScale) {
    Atoms atoms(
        {{0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}},
        {Species("Fe"), Species("Fe")}
    );

    SimpleRegularizationRules base_rules(0.001, 0.05, 0.1, 0.02);
    ScaledRegularizationRules scaled_rules(base_rules, /*force_scale=*/0.5);

    auto base_sig = base_rules.determine(atoms);
    auto scaled_sig = scaled_rules.determine(atoms);

    // With no forces on atoms, max force is 0.0, so scale is 1.0
    ASSERT_TRUE(base_sig.energy.has_value());
    ASSERT_TRUE(scaled_sig.energy.has_value());
    EXPECT_DOUBLE_EQ(*scaled_sig.energy, *base_sig.energy);

    ASSERT_TRUE(base_sig.forces.has_value());
    ASSERT_TRUE(scaled_sig.forces.has_value());
    EXPECT_DOUBLE_EQ((*scaled_sig.forces)[0].x, (*base_sig.forces)[0].x);

    ASSERT_TRUE(base_sig.virials.has_value());
    ASSERT_TRUE(scaled_sig.virials.has_value());
    EXPECT_DOUBLE_EQ(scaled_sig.virials->xx, base_sig.virials->xx);
}

TEST(TestScaledRegularizationRules, ScalesWithMaxForceMagnitude) {
    Atoms atoms(
        {{0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {2.0, 0.0, 0.0}},
        {Species("Fe"), Species("Fe"), Species("Fe")}
    );
    // Forces: atom 0 -> (3, 4, 0) magnitude 5.0; atom 1 -> (1, 0, 0) magnitude 1.0; atom 2 -> (0, 0, 0)
    atoms.setForces({{3.0, 4.0, 0.0}, {1.0, 0.0, 0.0}, {0.0, 0.0, 0.0}});

    SimpleRegularizationRules base_rules(0.002, 0.05, 0.1, 0.02);
    const Real force_scale = 0.4;
    // max_force = 5.0 => scale = 1.0 + 0.4 * 5.0 = 3.0
    ScaledRegularizationRules scaled_rules(base_rules, force_scale);

    auto base_sig = base_rules.determine(atoms);
    auto scaled_sig = scaled_rules.determine(atoms);

    const Real expected_scale = 3.0;

    ASSERT_TRUE(base_sig.energy.has_value());
    ASSERT_TRUE(scaled_sig.energy.has_value());
    EXPECT_NEAR(*scaled_sig.energy, *base_sig.energy * expected_scale, 1e-12);

    ASSERT_TRUE(base_sig.forces.has_value());
    ASSERT_TRUE(scaled_sig.forces.has_value());
    for (size_t i = 0; i < atoms.nAtoms(); ++i) {
        EXPECT_NEAR((*scaled_sig.forces)[i].x, (*base_sig.forces)[i].x * expected_scale, 1e-12);
        EXPECT_NEAR((*scaled_sig.forces)[i].y, (*base_sig.forces)[i].y * expected_scale, 1e-12);
        EXPECT_NEAR((*scaled_sig.forces)[i].z, (*base_sig.forces)[i].z * expected_scale, 1e-12);
    }

    ASSERT_TRUE(base_sig.virials.has_value());
    ASSERT_TRUE(scaled_sig.virials.has_value());
    EXPECT_NEAR(scaled_sig.virials->xx, base_sig.virials->xx * expected_scale, 1e-12);
    EXPECT_NEAR(scaled_sig.virials->xy, base_sig.virials->xy * expected_scale, 1e-12);
}

TEST(TestScaledRegularizationRules, RespectsMinScale) {
    Atoms atoms(
        {{0.0, 0.0, 0.0}},
        {Species("Fe")}
    );
    atoms.setForces({{0.0, 0.0, 0.0}});

    ScaledRegularizationRules scaled_rules(
        PerConfigTypeSigmas(0.001, 0.05, 0.1),
        /*force_scale=*/0.1,
        /*min_scale=*/2.5
    );

    auto sig = scaled_rules.determine(atoms);
    // Base energy sigma is 0.001, min_scale is 2.5 => scaled energy sigma = 0.0025
    ASSERT_TRUE(sig.energy.has_value());
    EXPECT_NEAR(*sig.energy, 0.0025, 1e-12);
}
