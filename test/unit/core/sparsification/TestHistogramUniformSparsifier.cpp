#include <gtest/gtest.h>
#include "jgap/core/sparsification/HistogramUniformSparsifier.hpp"

using namespace jgap;

TEST(TestHistogramUniformSparsifier, MinPointFiltersLowerBound) {
    // Descriptors with values ranging from 0.0 to 10.0
    std::vector<Descriptor<1>> descriptors;
    for (Real val = 0.0_r; val <= 10.0_r; val += 0.5_r) {
        descriptors.push_back({val});
    }

    // Sparsifier with a lower bound min_point = 2.0
    Descriptor<1> min_point{2.0_r};
    HistogramUniformSparsifier<1> sparsifier(
        42,
        5,
        std::nullopt,
        std::nullopt,
        min_point
    );

    auto sparse_points = sparsifier.selectSparsePoints(descriptors);
    ASSERT_FALSE(sparse_points.empty());
    EXPECT_EQ(sparse_points.size(), 5);

    // Verify all selected sparse points are >= min_point (2.0)
    for (const auto& pt: sparse_points) {
        EXPECT_GE(pt[0], 2.0_r);
    }
}
