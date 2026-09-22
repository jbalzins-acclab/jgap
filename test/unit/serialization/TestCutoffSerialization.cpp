#include <gtest/gtest.h>

#include "jgap/core/cutoff/CosCutoff.hpp"
#include "jgap/core/cutoff/PerriotPolynomialCutoff.hpp"
#include "jgap/serialization/SerializationRegistry.hpp"

using namespace jgap;

namespace {
    // Serialize to a temp file via the registry, then read it back, exercising the same code path the
    // apps use (create() stamps the format version, open() checks it).
    std::string tmpFile(const std::string& name) {
        return testing::TempDir() + "/" + name;
    }
}

TEST(TestCutoffSerialization, CosCutoff) {
    const auto file = tmpFile("cos_cutoff.h5");
    SerializationRegistry<CutoffFunction>::serialize(ValuePtr<CutoffFunction>(CosCutoff(4.5, 1.5)), file);

    auto rt = SerializationRegistry<CutoffFunction>::deserialize(file);
    auto restored = rt.as<CosCutoff>();
    ASSERT_NE(restored, nullptr);
    EXPECT_DOUBLE_EQ(restored->getCutoff(), 4.5);
    EXPECT_DOUBLE_EQ(restored->getCutoffTransitionWidth(), 1.5);
}

TEST(TestCutoffSerialization, PerriotPolynomialCutoff) {
    const auto file = tmpFile("perriot_cutoff.h5");
    SerializationRegistry<CutoffFunction>::serialize(
        ValuePtr<CutoffFunction>(PerriotPolynomialCutoff(6.0, 0.8)), file);

    auto rt = SerializationRegistry<CutoffFunction>::deserialize(file);
    auto restored = rt.as<PerriotPolynomialCutoff>();
    ASSERT_NE(restored, nullptr);
    EXPECT_DOUBLE_EQ(restored->getCutoff(), 6.0);
    EXPECT_DOUBLE_EQ(restored->getCutoffTransitionWidth(), 0.8);
}
