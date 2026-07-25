#include <gtest/gtest.h>

#include "jgap/core/kernels/SquaredExpKernel.hpp"
#include "jgap/experimental/kernels/CauchyKernel.hpp"
#include "jgap/experimental/kernels/WendlandKernel.hpp"
#include "jgap/serialization/SerializationRegistry.hpp"
#include "jgap/serialization/kernels/CauchyKernelSerialization.hpp"
#include "jgap/serialization/kernels/SquaredExpKernelSerialization.hpp"
#include "jgap/serialization/kernels/WendlandKernelSerialization.hpp"

using namespace jgap;

namespace {
    std::string tmpFile(const std::string& name) { return testing::TempDir() + "/" + name; }

    template<typename TKernel>
    TKernel roundTrip(const TKernel& kernel, const std::string& file) {
        ValuePtr<Kernel<TKernel::Dim>> kernel_ptr = kernel;
        SerializationRegistry<Kernel<TKernel::Dim>>::serialize(kernel_ptr, file);
        auto restored_ptr = SerializationRegistry<Kernel<TKernel::Dim>>::deserialize(file);
        auto* typed = restored_ptr.template as<TKernel>();
        if (!typed) {
            throw std::runtime_error("Failed to deserialize kernel to expected type");
        }
        return *typed;
    }

    template<typename TKernel>
    void expectKernelEq(const TKernel& a, const TKernel& b) {
        EXPECT_DOUBLE_EQ(a.getEnergyScale(), b.getEnergyScale());
        const auto la = a.getLengthScales();
        const auto lb = b.getLengthScales();
        for (size_t i = 0; i < TKernel::ExpDim; ++i) {
            EXPECT_DOUBLE_EQ(la[i], lb[i]);
        }
    }
}

TEST(TestKernelSerialization, SquaredExpKernel_1_0) {
    SquaredExpKernel<1, 0> kernel(2.5, {0.7});
    expectKernelEq(kernel, roundTrip(kernel, tmpFile("sqexp_1_0.h5")));
}

TEST(TestKernelSerialization, SquaredExpKernel_1_1) {
    SquaredExpKernel<1, 1> kernel(10.0, {1.3});
    expectKernelEq(kernel, roundTrip(kernel, tmpFile("sqexp_1_1.h5")));
}

TEST(TestKernelSerialization, SquaredExpKernel_3_1) {
    SquaredExpKernel<3, 1> kernel(1.5, {0.5, 1.0, 2.0});
    expectKernelEq(kernel, roundTrip(kernel, tmpFile("sqexp_3_1.h5")));
}

TEST(TestKernelSerialization, CauchyKernel_1_0) {
    CauchyKernel<1, 0> kernel(2.5, {0.7});
    expectKernelEq(kernel, roundTrip(kernel, tmpFile("cauchy_1_0.h5")));
}

TEST(TestKernelSerialization, CauchyKernel_1_1) {
    CauchyKernel<1, 1> kernel(10.0, {1.3});
    expectKernelEq(kernel, roundTrip(kernel, tmpFile("cauchy_1_1.h5")));
}

TEST(TestKernelSerialization, CauchyKernel_3_1) {
    CauchyKernel<3, 1> kernel(1.5, {0.5, 1.0, 2.0});
    expectKernelEq(kernel, roundTrip(kernel, tmpFile("cauchy_3_1.h5")));
}

TEST(TestKernelSerialization, WendlandKernel_1_0) {
    WendlandKernel<1, 0> kernel(2.5, {0.7});
    expectKernelEq(kernel, roundTrip(kernel, tmpFile("wendland_1_0.h5")));
}

TEST(TestKernelSerialization, WendlandKernel_1_1) {
    WendlandKernel<1, 1> kernel(2.0, {1.5});
    expectKernelEq(kernel, roundTrip(kernel, tmpFile("wendland_1_1.h5")));
}

TEST(TestKernelSerialization, WendlandKernel_3_1) {
    WendlandKernel<3, 1> kernel(1.5, {0.5, 1.0, 2.0});
    expectKernelEq(kernel, roundTrip(kernel, tmpFile("wendland_3_1.h5")));
}

TEST(TestKernelSerialization, WendlandKernelValues) {
    WendlandKernel<1, 0> kernel(3.0, {2.0}); // prefactor = 9.0, L = 2.0
    Descriptor<1> q1{0.0};
    Descriptor<1> q2{1.0}; // r = 1.0 / 2.0 = 0.5
    // val = 9.0 * (1 - 0.5)^4 * (1 + 4*0.5) = 9.0 * 0.0625 * 3.0 = 1.6875
    EXPECT_DOUBLE_EQ(kernel.value(q1, q2), 1.6875);

    Descriptor<1> q_far{3.0}; // r = 3.0 / 2.0 = 1.5 >= 1.0 => 0.0
    EXPECT_DOUBLE_EQ(kernel.value(q1, q_far), 0.0);
}
