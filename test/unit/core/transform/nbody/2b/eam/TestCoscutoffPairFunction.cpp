#include <cmath>
#include <gtest/gtest.h>

#include "jgap/core/atomic/geometry/Cluster2.hpp"
#include "jgap/core/transform/nbody/2b/eam/CoscutoffPairFunction.hpp"

using namespace jgap;

TEST(TestCoscutoffPairFunction, BelowRange) {
    CoscutoffPairFunction func(5.0, 2.0, 1.5); // cutoff, r_min, prefactor
    Cluster2 pair{.idx0 = 0, .idx1 = 1};
    pair.separation01.magnitude = 1.0;

    auto desc = func.evaluate(pair);
    EXPECT_DOUBLE_EQ(desc[0], 1.5);

    auto desc_deriv = func.evaluateAndDifferentiate(pair);
    EXPECT_DOUBLE_EQ(desc_deriv.value[0], 1.5);
    EXPECT_DOUBLE_EQ(desc_deriv.derivatives[0], 0.0);
}

TEST(TestCoscutoffPairFunction, AboveRange) {
    CoscutoffPairFunction func(5.0, 2.0, 1.5);
    Cluster2 pair{.idx0 = 0, .idx1 = 1};
    pair.separation01.magnitude = 6.0;

    auto desc = func.evaluate(pair);
    EXPECT_DOUBLE_EQ(desc[0], 0.0);

    auto desc_deriv = func.evaluateAndDifferentiate(pair);
    EXPECT_DOUBLE_EQ(desc_deriv.value[0], 0.0);
    EXPECT_DOUBLE_EQ(desc_deriv.derivatives[0], 0.0);
}

TEST(TestCoscutoffPairFunction, InsideRange) {
    CoscutoffPairFunction func(5.0, 2.0, 1.5); // r_min=2, cutoff=5, width=3
    Cluster2 pair{.idx0 = 0, .idx1 = 1};
    pair.separation01.magnitude = 3.5; // Midpoint, chi = 0.5, phase = pi/2

    // val = prefactor * 0.5 * (1 + cos(pi/2)) = 1.5 * 0.5 * (1 + 0) = 0.75
    auto desc = func.evaluate(pair);
    EXPECT_NEAR(desc[0], 0.75, 1e-14);

    // deriv = -prefactor * dchi_dr * 0.5 * pi * sin(pi/2)
    //       = -1.5 * (1/3) * 0.5 * pi * 1 = -0.25 * pi
    auto desc_deriv = func.evaluateAndDifferentiate(pair);
    EXPECT_NEAR(desc_deriv.value[0], 0.75, 1e-14);
    EXPECT_NEAR(desc_deriv.derivatives[0], -0.25 * M_PI, 1e-14);
}
