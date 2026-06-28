#include <gtest/gtest.h>
#include "core/atomic/io/XYZData.hpp"
#include "core/atomic/neighbours/NeighbourList.hpp"
#include "core/atomic/Atoms.hpp"

using namespace jgap;

TEST(TestNeighbourList, Sample1CorrectNumberOfNeighbours) {
    auto atoms = Atoms::readAtoms("test/resources/xyz-samples/iter-3-3-test.xyz");
    auto neighbour_lists = NeighbourList::form(atoms, 5.0);

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

TEST(TestNeighbourList, Sample2CorrectNumberOfNeighbours) {
    auto atoms = Atoms::readAtoms("test/resources/xyz-samples/iter-3-3-train.xyz");
    auto neighbour_lists = NeighbourList::form(atoms, 5.0);

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

TEST(TestNeighbourList, FindAllClusters) {
    // A square of 4 atoms. Side length 10. Diagonal is ~14.14.
    // Cutoff is 11, so adjacent atoms are neighbors, but diagonal ones are not.
    Atoms atoms(
        { {0,0,0}, {10,0,0}, {0,10,0}, {10,10,0} },
        { Species("Fe"), Species("Ni"), Species("Fe"), Species("Ni") }
    );
    auto nl = NeighbourList(atoms, 11.0);

    // --- Test findAllClusters<2> ---
    // Each atom has 2 neighbors. Total pairs = 4 atoms * 2 neighbors = 8.

    // Fe-Fe pairs: Fe(0)-Fe(2) and Fe(2)-Fe(0)
    auto fe_fe_pairs = nl.findAllClusters<ValueOnly>(SpeciesSet<2, NodeSymmetric>{"Fe", "Fe"});
    EXPECT_EQ(fe_fe_pairs.size(), 2);

    // Ni-Ni pairs: Ni(1)-Ni(3) and Ni(3)-Ni(1)
    auto ni_ni_pairs = nl.findAllClusters<ValueOnly>(SpeciesSet<2, NodeSymmetric>{"Ni", "Ni"});
    EXPECT_EQ(ni_ni_pairs.size(), 2);

    // Fe-Ni pairs: Fe(0)-Ni(1) and Fe(2)-Ni(3)
    auto fe_ni_pairs = nl.findAllClusters<ValueOnly>(SpeciesSet<2, NodeSymmetric>{"Fe", "Ni"});
    EXPECT_EQ(fe_ni_pairs.size(), 2);

    // Ni-Fe pairs: Ni(1)-Fe(0) and Ni(3)-Fe(2)
    auto ni_fe_pairs = nl.findAllClusters(SpeciesSet<2, NodeSymmetric>{"Ni", "Fe"});
    EXPECT_EQ(ni_fe_pairs.size(), 2);

    auto fe_ni_pairs_symm = nl.findAllClusters(SpeciesSet<2, FullSymmetry>{"Fe", "Ni"});
    EXPECT_EQ(fe_ni_pairs_symm.size(), 2);

    auto ni_fe_pairs_symm = nl.findAllClusters(SpeciesSet<2, FullSymmetry>{"Ni", "Fe"});
    EXPECT_EQ(ni_fe_pairs_symm.size(), 2);

    // --- Test findAllClusters<3> ---
    // Triplets are formed from a central atom and two of its neighbors.

    // {"Fe", "Fe", "Ni"}:
    // Centered on Fe(0), neighbors are Fe(2) and Ni(1). Forms triplet (0,2,1).
    // Centered on Fe(2), neighbors are Fe(0) and Ni(3). Forms triplet (2,0,3).
    auto fe_fe_ni_triplets = nl.findAllClusters(SpeciesSet<3, NodeSymmetric>{"Fe", "Fe", "Ni"});
    EXPECT_EQ(fe_fe_ni_triplets.size(), 2);

    // {"Ni", "Ni", "Fe"}:
    // Centered on Ni(1), neighbors are Ni(3) and Fe(0). Forms triplet (1,3,0).
    // Centered on Ni(3), neighbors are Ni(1) and Fe(2). Forms triplet (3,1,2).
    auto ni_ni_fe_triplets = nl.findAllClusters(SpeciesSet<3, NodeSymmetric>{"Ni", "Ni", "Fe"});
    EXPECT_EQ(ni_ni_fe_triplets.size(), 2);

    // Other combinations like {"Fe", "Ni", "Ni"} are impossible as no Fe has two Ni neighbors.
    auto fe_ni_ni_triplets = nl.findAllClusters(SpeciesSet<3, NodeSymmetric>{"Fe", "Ni", "Ni"});
    EXPECT_EQ(fe_ni_ni_triplets.size(), 0);
}