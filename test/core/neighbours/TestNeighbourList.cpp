#include <fstream>
#include <gtest/gtest.h>

#include <core/neighbours/NeighbourFinder.hpp>

#include "utils/Utils.hpp"

TEST(NeighbourFinderTest, Sample1CorrectNumberOfNeighbours) {

    // TODO: fix path
    auto structs = jgap::readXyz("test/resources/xyz-samples/iter-3-3-test.xyz");

    jgap::NeighbourFinder::findNeighbours(structs, 5);

    size_t total = 0;
    for (const auto &structure: structs) {
        for (auto &atomData : structure.atoms) {
            total += atomData.neighbours.value().size();
        }
    }
    EXPECT_EQ(total, 862046);
}

TEST(NeighbourFinderTest, Sample2CorrectNumberOfNeighbours) {

    // TODO: fix path
    auto structs = jgap::readXyz("test/resources/xyz-samples/iter-3-3-train.xyz");

    jgap::NeighbourFinder::findNeighbours(structs, 5);

    size_t total = 0;
    for (const auto &structure: structs) {
        for (auto &atomData : structure.atoms) {
            total += atomData.neighbours.value().size();
        }
    }
    EXPECT_EQ(total, 8815680);
}