#include <cstdio>
#include <gtest/gtest.h>
#include "jgap/experimental/cutoff/WendlandFunction.hpp"
#include "jgap/experimental/transform/nbody/2b/CoordinationTransformation.hpp"
#include "jgap/serialization/SerializationRegistry.hpp"

using namespace jgap;

TEST(TestCoordinationSerialization, WendlandFunctionSerialization) {
    auto wendland = ValuePtr<CutoffFunction>(WendlandFunction(1.5, 3.2));

    std::string temp_filename = "test_wendland_serialization.h5";
    SerializationRegistry<CutoffFunction>::serialize(wendland, temp_filename);

    auto deserialized = SerializationRegistry<CutoffFunction>::deserialize(temp_filename);

    ASSERT_TRUE(deserialized);
    auto typed_deserialized = dynamic_cast<WendlandFunction*>(deserialized.get());
    ASSERT_NE(typed_deserialized, nullptr);

    EXPECT_DOUBLE_EQ(typed_deserialized->getRMin(), 1.5);
    EXPECT_DOUBLE_EQ(typed_deserialized->getRMax(), 3.2);
    EXPECT_DOUBLE_EQ(typed_deserialized->getCutoff(), 3.2);

    std::remove(temp_filename.c_str());
}

TEST(TestCoordinationsSerialization, CoordinationsTransformationSerialization) {
    std::array<std::pair<double, double>, 2> ranges = {{{2.4, 2.6}, {2.75, 2.95}}};
    auto transform = ValuePtr<TwoBodyTransformation<2>>(CoordinationTransformation<2>(ranges));

    std::string temp_filename = "test_coordinations_transform_serialization.h5";
    SerializationRegistry<TwoBodyTransformation<2>>::serialize(transform, temp_filename);

    auto deserialized = SerializationRegistry<TwoBodyTransformation<2>>::deserialize(temp_filename);

    ASSERT_TRUE(deserialized);
    auto typed_deserialized = dynamic_cast<CoordinationTransformation<2>*>(deserialized.get());
    ASSERT_NE(typed_deserialized, nullptr);

    auto deserialized_ranges = typed_deserialized->getRanges();
    EXPECT_DOUBLE_EQ(deserialized_ranges[0].first, 2.4);
    EXPECT_DOUBLE_EQ(deserialized_ranges[0].second, 2.6);
    EXPECT_DOUBLE_EQ(deserialized_ranges[1].first, 2.75);
    EXPECT_DOUBLE_EQ(deserialized_ranges[1].second, 2.95);

    std::remove(temp_filename.c_str());
}
