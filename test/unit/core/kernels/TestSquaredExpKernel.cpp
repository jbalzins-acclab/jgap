#include <cmath>
#include <gtest/gtest.h>
#include "jgap/core/kernels/SquaredExpKernel.hpp"

using namespace jgap;

TEST(TestSquaredExpKernel, Kernel_1_0) {
    // E=1, C=0. Pure exponential kernel.
    SquaredExpKernel<1, 0> kernel(2.0, {0.5}); // energy_scale, length_scales
    std::array<Real, 1> q1 = {2.0};
    std::array<Real, 1> q2 = {3.0};

    // prefactor = 2^2 = 4.0
    // inv_l_sq[0] = 1 / 0.5^2 = 4.0
    // exp_arg = (2-3)^2 * 4.0 = 4.0
    // value = 4.0 * exp(-0.5 * 4.0) = 4.0 * exp(-2.0)
    Real expected_value = 4.0 * std::exp(-2.0);
    EXPECT_NEAR(kernel.value(q1, q2), expected_value, 1e-9);

    // grad[0] = value * (q1[0]-q2[0]) * inv_l_sq[0]
    //         = expected_value * (2.0-3.0) * 4.0 = -4.0 * expected_value
    Real expected_grad0 = -4.0 * expected_value;
    auto [val, grad] = kernel.valueAndGradient(q1, q2);
    EXPECT_NEAR(val, expected_value, 1e-9);
    EXPECT_NEAR(grad[0], expected_grad0, 1e-9);
}

TEST(TestSquaredExpKernel, Kernel_1_1) {
    // E=1, C=1.
    SquaredExpKernel<1, 1> kernel(2.0, {0.5});
    std::array<Real, 2> q1 = {2.0, 0.5};
    std::array<Real, 2> q2 = {3.0, 0.8};

    // exp_val = 4.0 * exp(-2.0)
    // value = exp_val * (q1[1]*q2[1]) = exp_val * (0.5*0.8) = exp_val * 0.4
    Real expected_value = 4.0 * std::exp(-2.0) * 0.4;
    EXPECT_NEAR(kernel.value(q1, q2), expected_value, 1e-9);

    // grad[0] = value * (q1[0]-q2[0]) * inv_l_sq[0] = expected_value * -1 * 4
    Real expected_grad0 = -4.0 * expected_value;
    // grad[1] = value / q2[1] = expected_value / 0.8
    Real expected_grad1 = expected_value / 0.8;

    auto [val, grad] = kernel.valueAndGradient(q1, q2);
    EXPECT_NEAR(val, expected_value, 1e-9);
    EXPECT_NEAR(grad[0], expected_grad0, 1e-9);
    EXPECT_NEAR(grad[1], expected_grad1, 1e-9);
}

TEST(TestSquaredExpKernel, Kernel_2_1) {
    // E=2, C=1.
    SquaredExpKernel<2, 1> kernel(2.0, {0.5, 1.0});
    std::array<Real, 3> q1 = {2.0, 1.0, 0.5};
    std::array<Real, 3> q2 = {3.0, 1.5, 0.8};

    // exp_arg = (2-3)^2*4 + (1-1.5)^2*1 = 4.0 + 0.25 = 4.25
    // exp_val = 4.0 * exp(-0.5 * 4.25) = 4.0 * exp(-2.125)
    // value = exp_val * (0.5 * 0.8) = exp_val * 0.4
    Real expected_value = 4.0 * std::exp(-2.125) * 0.4;
    EXPECT_NEAR(kernel.value(q1, q2), expected_value, 1e-9);

    // grad[0] = value * (2-3) * 4 = -4 * value
    Real expected_grad0 = -4.0 * expected_value;
    // grad[1] = value * (1-1.5) * 1 = -0.5 * value
    Real expected_grad1 = -0.5 * expected_value;
    // grad[2] = value / 0.8
    Real expected_grad2 = expected_value / 0.8;

    auto [val, grad] = kernel.valueAndGradient(q1, q2);
    EXPECT_NEAR(val, expected_value, 1e-9);
    EXPECT_NEAR(grad[0], expected_grad0, 1e-9);
    EXPECT_NEAR(grad[1], expected_grad1, 1e-9);
    EXPECT_NEAR(grad[2], expected_grad2, 1e-9);
}

TEST(TestSquaredExpKernel, Kernel_3_1) {
    // E=3, C=1.
    SquaredExpKernel<3, 1> kernel(1.0, {1.0, 1.0, 1.0});
    std::array<Real, 4> q1 = {1.0, 2.0, 3.0, 0.5};
    std::array<Real, 4> q2 = {1.5, 2.5, 3.5, 0.8};

    // exp_arg = (0.5^2 + 0.5^2 + 0.5^2) * 1 = 3 * 0.25 = 0.75
    // exp_val = 1.0 * exp(-0.5 * 0.75) = exp(-0.375)
    // value = exp_val * (0.5 * 0.8) = 0.4 * exp(-0.375)
    Real expected_value = 0.4 * std::exp(-0.375);
    EXPECT_NEAR(kernel.value(q1, q2), expected_value, 1e-9);

    // grad[0] = value * (1-1.5) * 1 = -0.5 * value
    Real expected_grad0 = -0.5 * expected_value;
    // grad[1] = value * (2-2.5) * 1 = -0.5 * value
    Real expected_grad1 = -0.5 * expected_value;
    // grad[2] = value * (3-3.5) * 1 = -0.5 * value
    Real expected_grad2 = -0.5 * expected_value;
    // grad[3] = value / 0.8
    Real expected_grad3 = expected_value / 0.8;

    auto [val, grad] = kernel.valueAndGradient(q1, q2);
    EXPECT_NEAR(val, expected_value, 1e-9);
    EXPECT_NEAR(grad[0], expected_grad0, 1e-9);
    EXPECT_NEAR(grad[1], expected_grad1, 1e-9);
    EXPECT_NEAR(grad[2], expected_grad2, 1e-9);
    EXPECT_NEAR(grad[3], expected_grad3, 1e-9);
}
