#include <cmath>
#include <gtest/gtest.h>

#include "jgap/core/atomic/geometry/Cluster2.hpp"
#include "jgap/core/transform/nbody/2b/eam/FSGenPairFunction.hpp"

using namespace jgap;

TEST(TestFSGenPairFunction, AtOrigin) {
    FSGenPairFunction func(5.0, 3.0, 1.5); // cutoff, degree, prefactor
    Cluster2 pair{.atom_indexes = {0, 1}, .r01 = 0.0};

    auto desc = func.evaluate(pair);
    EXPECT_DOUBLE_EQ(desc[0], 1.5); // prefactor * (1 - 0)^3

    auto desc_deriv = func.evaluateAndDifferentiate(pair);
    EXPECT_DOUBLE_EQ(desc_deriv.value[0], 1.5);

    // deriv = -prefactor * (1 - 0)^2 * degree / cutoff = -1.5 * 1 * 3 / 5 = -0.9
    EXPECT_DOUBLE_EQ(desc_deriv.derivatives[0], -0.9);
}

TEST(TestFSGenPairFunction, AboveCutoff) {
    FSGenPairFunction func(5.0, 3.0, 1.5);
    Cluster2 pair{.atom_indexes = {0, 1}, .r01 = 6.0};

    auto desc = func.evaluate(pair);
    EXPECT_DOUBLE_EQ(desc[0], 0.0);

    auto desc_deriv = func.evaluateAndDifferentiate(pair);
    EXPECT_DOUBLE_EQ(desc_deriv.value[0], 0.0);
    EXPECT_DOUBLE_EQ(desc_deriv.derivatives[0], 0.0);
}

TEST(TestFSGenPairFunction, InsideRange) {
    FSGenPairFunction func(5.0, 3.0, 1.5);
    Cluster2 pair{.atom_indexes = {0, 1}, .r01 = 2.5}; // midpoint, distance * cutoff_inverse = 0.5

    // val = prefactor * (1 - 0.5)^3 = 1.5 * 0.125 = 0.1875
    auto desc = func.evaluate(pair);
    EXPECT_DOUBLE_EQ(desc[0], 0.1875);

    // deriv = -prefactor * (1 - 0.5)^2 * degree * cutoff_inverse
    //       = -1.5 * 0.25 * 3 * 0.2 = -0.225
    auto desc_deriv = func.evaluateAndDifferentiate(pair);
    EXPECT_DOUBLE_EQ(desc_deriv.value[0], 0.1875);
    EXPECT_DOUBLE_EQ(desc_deriv.derivatives[0], -0.225);
}
