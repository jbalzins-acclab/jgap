#include <gtest/gtest.h>

#include "jgap/core/atomic/geometry/Cluster2.hpp"
#include "jgap/core/transform/nbody/2b/eam/PolycutoffPairFunction.hpp"

using namespace jgap;

TEST(TestPolycutoffPairFunction, BelowRange) {
    PolycutoffPairFunction func(5.0, 2.0, 1.5); // cutoff, r_min, prefactor
    Cluster2 pair{.idx0 = 0, .idx1 = 1};
    pair.separation01.magnitude = 1.0;
    pair.separation01.direction = Vector3{1.0, 0.0, 0.0};

    auto desc = func.evaluate(pair);
    EXPECT_DOUBLE_EQ(desc[0], 1.5);

    auto desc_deriv = func.evaluateAndDifferentiate(pair);
    EXPECT_DOUBLE_EQ(desc_deriv.value[0], 1.5);
    EXPECT_DOUBLE_EQ(desc_deriv.grad_r1[0].x, 0.0);
}

TEST(TestPolycutoffPairFunction, AboveRange) {
    PolycutoffPairFunction func(5.0, 2.0, 1.5);
    Cluster2 pair{.idx0 = 0, .idx1 = 1};
    pair.separation01.magnitude = 6.0;
    pair.separation01.direction = Vector3{1.0, 0.0, 0.0};

    auto desc = func.evaluate(pair);
    EXPECT_DOUBLE_EQ(desc[0], 0.0);

    auto desc_deriv = func.evaluateAndDifferentiate(pair);
    EXPECT_DOUBLE_EQ(desc_deriv.value[0], 0.0);
    EXPECT_DOUBLE_EQ(desc_deriv.grad_r1[0].x, 0.0);
}

TEST(TestPolycutoffPairFunction, InsideRange) {
    PolycutoffPairFunction func(5.0, 2.0, 1.5); // r_min=2, cutoff=5, width=3
    Cluster2 pair{.idx0 = 0, .idx1 = 1};
    pair.separation01.magnitude = 3.5; // midpoint, chi = 0.5
    pair.separation01.direction = Vector3{1.0, 0.0, 0.0};

    // val = prefactor * (1.0 - chi^3 * (6*chi^2 - 15*chi + 10))
    // chi = 0.5 -> poly part is 0.5 (as proven in Perriot test)
    // val = 1.5 * 0.5 = 0.75
    auto desc = func.evaluate(pair);
    EXPECT_DOUBLE_EQ(desc[0], 0.75);

    // deriv = prefactor * dchi_dr * chi^2 * (-30*chi^2 + 60*chi - 30)
    // dchi_dr = 1/3
    // chi = 0.5 -> (-30*0.25 + 60*0.5 - 30) = -7.5 + 30 - 30 = -7.5
    // deriv = 1.5 * (1/3) * 0.25 * (-7.5) = 0.5 * 0.25 * (-7.5) = 0.125 * -7.5 = -0.9375
    auto desc_deriv = func.evaluateAndDifferentiate(pair);
    EXPECT_DOUBLE_EQ(desc_deriv.value[0], 0.75);
    EXPECT_DOUBLE_EQ(desc_deriv.grad_r1[0].x, -0.9375);
}
