#include <gtest/gtest.h>
#include "core/cutoff/CosCutoff.hpp"
#include <cmath>

using namespace jgap;

TEST(TestCosCutoff, BelowCutoffRange) {
    CosCutoff cutoff(5.0, 3.0); // r_min = 2.0

    // r < r_min
    EXPECT_DOUBLE_EQ(cutoff.evaluate(1.0), 1.0);

    auto [val, deriv] = cutoff.evaluateAndDifferentiate(1.0);
    EXPECT_DOUBLE_EQ(val, 1.0);
    EXPECT_DOUBLE_EQ(deriv, 0.0);
}

TEST(TestCosCutoff, AboveCutoffRange) {
    CosCutoff cutoff(5.0, 3.0); // r_min = 2.0

    // r > cutoff
    EXPECT_DOUBLE_EQ(cutoff.evaluate(6.0), 0.0);

    auto [val, deriv] = cutoff.evaluateAndDifferentiate(6.0);
    EXPECT_DOUBLE_EQ(val, 0.0);
    EXPECT_DOUBLE_EQ(deriv, 0.0);
}

TEST(TestCosCutoff, InsideCutoffRange) {
    CosCutoff cutoff(5.0, 3.0); // r_min = 2.0, cutoff = 5.0, width = 3.0

    // r = 3.5 (midpoint) -> phase = (3.5 - 2.0) * pi / 3.0 = 1.5 * pi / 3.0 = 0.5 * pi
    // val = 0.5 * (cos(pi/2) + 1.0) = 0.5 * (0 + 1) = 0.5
    // deriv = -0.5 * (pi / 3.0) * sin(pi/2) = -0.5 * pi / 3.0 = -pi / 6.0

    EXPECT_NEAR(cutoff.evaluate(3.5), 0.5, 1e-14);

    auto [val, deriv] = cutoff.evaluateAndDifferentiate(3.5);
    EXPECT_NEAR(val, 0.5, 1e-14);
    EXPECT_NEAR(deriv, -M_PI / 6.0, 1e-14);

    // r = 2.0 (boundary) -> phase = 0
    // val = 1.0
    // deriv = 0.0
    auto [val_min, deriv_min] = cutoff.evaluateAndDifferentiate(2.0);
    EXPECT_DOUBLE_EQ(val_min, 1.0);
    EXPECT_DOUBLE_EQ(deriv_min, 0.0);

    // r = 5.0 (boundary) -> phase = pi
    // val = 0.0
    // deriv = 0.0
    auto [val_max, deriv_max] = cutoff.evaluateAndDifferentiate(5.0);
    EXPECT_NEAR(val_max, 0.0, 1e-14);
    EXPECT_NEAR(deriv_max, 0.0, 1e-14);
}
