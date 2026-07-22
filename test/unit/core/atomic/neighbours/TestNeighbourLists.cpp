#include <gtest/gtest.h>
#include "core/atomic/io/XYZData.hpp"
#include "core/atomic/neighbours/NeighbourLists.hpp"
#include "core/atomic/Atoms.hpp"

using namespace jgap;

TEST(TestNeighbourLists, Sample1CorrectNumberOfNeighbours) {
    auto atoms = Atoms::readAtoms("test/resources/xyz-samples/iter-3-3-test.xyz");
    auto neighbour_lists = NeighbourLists::form(atoms, 5.0);

    size_t total_neighbours = 0;
    for (const auto& nl : neighbour_lists) {
        for (const auto& atom_neighbours : nl.neighbours_per_atom) {
            for (const auto& [species, neighbours] : atom_neighbours) {
                total_neighbours += neighbours.size();
            }
        }
    }
    EXPECT_EQ(total_neighbours, 862046);
}

TEST(TestNeighbourLists, Sample2CorrectNumberOfNeighbours) {
    auto atoms = Atoms::readAtoms("test/resources/xyz-samples/iter-3-3-train.xyz");
    auto neighbour_lists = NeighbourLists::form(atoms, 5.0);

    size_t total_neighbours = 0;
    for (const auto& nl : neighbour_lists) {
        for (const auto& atom_neighbours : nl.neighbours_per_atom) {
            for (const auto& [species, neighbours] : atom_neighbours) {
                total_neighbours += neighbours.size();
            }
        }
    }
    EXPECT_EQ(total_neighbours, 8815680);
}
