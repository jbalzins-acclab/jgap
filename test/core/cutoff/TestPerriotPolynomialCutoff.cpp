#include <gtest/gtest.h>
#include "core/cutoff/PerriotPolynomialCutoff.hpp"

using namespace jgap;

TEST(TestPerriotPolynomialCutoff, BelowCutoffRange) {
    PerriotPolynomialCutoff cutoff(2.0, 5.0);

    // r < r_min
    EXPECT_DOUBLE_EQ(cutoff.evaluate(1.0), 1.0);

    auto [val, deriv] = cutoff.evaluateAndDifferentiate(1.0);
    EXPECT_DOUBLE_EQ(val, 1.0);
    EXPECT_DOUBLE_EQ(deriv, 0.0);
}

TEST(TestPerriotPolynomialCutoff, AboveCutoffRange) {
    PerriotPolynomialCutoff cutoff(2.0, 5.0);

    // r > cutoff
    EXPECT_DOUBLE_EQ(cutoff.evaluate(6.0), 0.0);

    auto [val, deriv] = cutoff.evaluateAndDifferentiate(6.0);
    EXPECT_DOUBLE_EQ(val, 0.0);
    EXPECT_DOUBLE_EQ(deriv, 0.0);
}

TEST(TestPerriotPolynomialCutoff, InsideCutoffRange) {
    PerriotPolynomialCutoff cutoff(2.0, 5.0); // r_min = 2.0, cutoff = 5.0, width = 3.0

    // r = 3.5 (midpoint) -> chi = (3.5 - 2.0) / 3.0 = 0.5
    // val = 1.0 - 0.5^3 * (6(0.25) - 15(0.5) + 10) = 1.0 - 0.125 * (1.5 - 7.5 + 10) = 1.0 - 0.125 * 4 = 1.0 - 0.5 = 0.5
    // deriv = -30 * 0.25 * (-0.5)^2 / 3.0 = -30 * 0.25 * 0.25 / 3.0 = -30 * 0.0625 / 3.0 = -1.875 / 3.0 = -0.625

    EXPECT_DOUBLE_EQ(cutoff.evaluate(3.5), 0.5);

    auto [val, deriv] = cutoff.evaluateAndDifferentiate(3.5);
    EXPECT_DOUBLE_EQ(val, 0.5);
    EXPECT_DOUBLE_EQ(deriv, -0.625);

    // r = 2.0 (boundary) -> chi = 0
    // val = 1.0
    // deriv = 0.0
    auto [val_min, deriv_min] = cutoff.evaluateAndDifferentiate(2.0);
    EXPECT_DOUBLE_EQ(val_min, 1.0);
    EXPECT_DOUBLE_EQ(deriv_min, 0.0);

    // r = 5.0 (boundary) -> chi = 1
    // val = 0.0
    // deriv = 0.0
    auto [val_max, deriv_max] = cutoff.evaluateAndDifferentiate(5.0);
    EXPECT_DOUBLE_EQ(val_max, 0.0);
    EXPECT_DOUBLE_EQ(deriv_max, 0.0);
}
