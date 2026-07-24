#include <gtest/gtest.h>

#include "jgap/core/kernels/SquaredExpKernel.hpp"
#include "jgap/serialization/SerializationNode.hpp"
#include "jgap/serialization/kernels/SquaredExpKernelSerialization.hpp"

using namespace jgap;

namespace {
    std::string tmpFile(const std::string& name) { return testing::TempDir() + "/" + name; }

    // SquaredExpKernel isn't registry-managed (it's a value type keyed on its dimensions), so it is
    // round-tripped through a node directly via SquaredExpKernelSerialization<E, C>.
    template<size_t E, size_t C>
    SquaredExpKernel<E, C> roundTrip(const SquaredExpKernel<E, C>& kernel, const std::string& file) {
        {
            SerializationNode node = SerializationNode::create(file);
            SquaredExpKernelSerialization<E, C>::serialize(kernel, node);
        }
        SerializationNode node = SerializationNode::open(file);
        return SquaredExpKernelSerialization<E, C>::deserialize(node);
    }

    template<size_t E, size_t C>
    void expectKernelEq(const SquaredExpKernel<E, C>& a, const SquaredExpKernel<E, C>& b) {
        EXPECT_DOUBLE_EQ(a.getEnergyScale(), b.getEnergyScale());
        const auto la = a.getLengthScales();
        const auto lb = b.getLengthScales();
        for (size_t i = 0; i < E; ++i) {
            EXPECT_DOUBLE_EQ(la[i], lb[i]);
        }
    }
}

TEST(TestKernelSerialization, Kernel_1_0) {
    SquaredExpKernel<1, 0> kernel(2.5, {0.7});
    expectKernelEq(kernel, roundTrip(kernel, tmpFile("kernel_1_0.h5")));
}

TEST(TestKernelSerialization, Kernel_1_1) {
    SquaredExpKernel<1, 1> kernel(10.0, {1.3});
    expectKernelEq(kernel, roundTrip(kernel, tmpFile("kernel_1_1.h5")));
}

TEST(TestKernelSerialization, Kernel_3_1) {
    SquaredExpKernel<3, 1> kernel(1.5, {0.5, 1.0, 2.0});
    expectKernelEq(kernel, roundTrip(kernel, tmpFile("kernel_3_1.h5")));
}
