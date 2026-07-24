#include <gtest/gtest.h>
#include "jgap/core/atomic/Atoms.hpp"
#include "jgap/core/atomic/iteration/AtomicCluster2Expansion.hpp"
#include "jgap/core/atomic/iteration/AtomicCluster3Expansion.hpp"
#include "jgap/core/atomic/iteration/Cluster2Expansion.hpp"
#include "jgap/core/atomic/neighbours/NeighbourLists.hpp"
#include "jgap/experimental/atomic/iteration/Cluster3Expansion.hpp"

using namespace jgap;

TEST(TestClusterExpansions, PBC1AtomCell) {
    Atoms atoms({{0.0, 0.0, 0.0}}, {Species("Fe")}, Lattice({2.0, 0.0, 0.0}, {0.0, 2.0, 0.0}, {0.0, 0.0, 2.0}),
                {true, false, false});
    // cutoff 2.5 means it will see its own periodic images at +2.0 and -2.0.
    auto nl = NeighbourLists(atoms, 2.5);

    // Cluster2Expansion gives 1 unique pair (tie-broken by direction.isPositive()).
    Cluster2Expansion c2e(Species2Sorted("Fe,Fe"));
    EXPECT_EQ(c2e.find(nl, CalculationType::ValueOnly).clusters.size(), 1);

    // AtomicCluster2Expansion gives 2 atom-centric pairs (both directions).
    AtomicCluster2Expansion ac2e(Species2Atomic("Fe|Fe"));
    EXPECT_EQ(ac2e.find(0, nl, CalculationType::ValueOnly).clusters.size(), 2);
}

TEST(TestClusterExpansions, PBC1AtomCellCluster3) {
    Atoms atoms({{0.0, 0.0, 0.0}}, {Species("Fe")}, Lattice({2.0, 0.0, 0.0}, {0.0, 2.0, 0.0}, {0.0, 0.0, 2.0}),
                {true, true, false});

    auto nl = NeighbourLists(atoms, 2.5);

    Cluster3Expansion c3e(Species3Sorted("Fe,Fe,Fe"));
    // jk - cutoff
    EXPECT_EQ(c3e.find(nl, CalculationType::ValueOnly).clusters.size(), 0);

    // AtomicCluster3Expansion gives 4 atom-centric triplets
    AtomicCluster3Expansion ac3e(Species3AtomicSorted("Fe|Fe,Fe"));
    EXPECT_EQ(ac3e.find(0, nl, CalculationType::ValueOnly).clusters.size(), 6);
}

TEST(TestClusterExpansions, Cluster2Expansion) {
    // A square of 4 atoms. Side length 10. Diagonal is ~14.14.
    // Cutoff is 11, so adjacent atoms are neighbors, but diagonal ones are not.
    Atoms atoms({{0, 0, 0}, {10, 0, 0}, {0, 10, 0}, {10, 10, 0}},
                {Species("Fe"), Species("Ni"), Species("Fe"), Species("Ni")});
    auto nl = NeighbourLists(atoms, 11.0);

    // Each atom has 2 neighbors. Total pairs = 4 atoms * 2 neighbors / 2 = 4 unique pairs.

    // Fe-Fe pair: Fe(0)-Fe(2)
    Cluster2Expansion fe_fe(Species2Sorted("Fe,Fe"));
    auto fe_fe_pairs = fe_fe.find(nl, CalculationType::ValueOnly).clusters;
    EXPECT_EQ(fe_fe_pairs.size(), 1);

    // Ni-Ni pair: Ni(1)-Ni(3)
    Cluster2Expansion ni_ni(Species2Sorted("Ni,Ni"));
    auto ni_ni_pairs = ni_ni.find(nl, CalculationType::ValueOnly).clusters;
    EXPECT_EQ(ni_ni_pairs.size(), 1);

    // Fe-Ni pairs: Fe(0)-Ni(1) and Fe(2)-Ni(3)
    Cluster2Expansion fe_ni(Species2Sorted("Fe,Ni"));
    auto fe_ni_pairs = fe_ni.find(nl, CalculationType::ValueOnly).clusters;
    EXPECT_EQ(fe_ni_pairs.size(), 2);
}

TEST(TestClusterExpansions, AtomicCluster2Expansion) {
    // A square of 4 atoms. Side length 10. Diagonal is ~14.14.
    // Cutoff is 11, so adjacent atoms are neighbors, but diagonal ones are not.
    Atoms atoms({{0, 0, 0}, {10, 0, 0}, {0, 10, 0}, {10, 10, 0}},
                {Species("Fe"), Species("Ni"), Species("Fe"), Species("Ni")});
    auto nl = NeighbourLists(atoms, 11.0);

    // Atom-centric pairs are formed from a central atom and its neighbors.

    // {"Fe", "Fe"}:
    // Centered on Fe(0), neighbor is Fe(2).
    AtomicCluster2Expansion fe_fe(Species2Atomic("Fe|Fe"));
    auto fe_fe_pairs = fe_fe.find(0, nl, CalculationType::ValueOnly).clusters;
    EXPECT_EQ(fe_fe_pairs.size(), 1);

    // {"Fe", "Ni"}:
    // Centered on Fe(0), neighbor is Ni(1).
    AtomicCluster2Expansion fe_ni(Species2Atomic("Fe|Ni"));
    auto fe_ni_pairs = fe_ni.find(0, nl, CalculationType::ValueOnly).clusters;
    EXPECT_EQ(fe_ni_pairs.size(), 1);

    // {"Ni", "Fe"}:
    // Centered on Ni(1), neighbor is Fe(0).
    AtomicCluster2Expansion ni_fe(Species2Atomic("Ni|Fe"));
    auto ni_fe_pairs = ni_fe.find(1, nl, CalculationType::ValueOnly).clusters;
    EXPECT_EQ(ni_fe_pairs.size(), 1);
}

TEST(TestClusterExpansions, AtomicCluster3Expansion) {
    // Same square of 4 atoms.
    Atoms atoms({{0, 0, 0}, {10, 0, 0}, {0, 10, 0}, {10, 10, 0}},
                {Species("Fe"), Species("Ni"), Species("Fe"), Species("Ni")});
    auto nl = NeighbourLists(atoms, 11.0);

    // Triplets are formed from a central atom and two of its neighbors.

    // {"Fe", "Fe", "Ni"}:
    // Centered on Fe(0), neighbors are Fe(2) and Ni(1). Forms triplet (0,2,1).
    // Centered on Fe(2), neighbors are Fe(0) and Ni(3). Forms triplet (2,0,3).
    AtomicCluster3Expansion fe_fe_ni(Species3AtomicSorted("Fe|Fe,Ni"));
    auto fe_fe_ni_triplets = fe_fe_ni.findAll(nl, CalculationType::ValueOnly).clusters;
    EXPECT_EQ(fe_fe_ni_triplets.size(), 2);

    // {"Ni", "Ni", "Fe"}:
    // Centered on Ni(1), neighbors are Ni(3) and Fe(0). Forms triplet (1,3,0).
    // Centered on Ni(3), neighbors are Ni(1) and Fe(2). Forms triplet (3,1,2).
    AtomicCluster3Expansion ni_ni_fe(Species3AtomicSorted("Ni|Fe,Ni"));
    auto ni_ni_fe_triplets = ni_ni_fe.findAll(nl, CalculationType::ValueOnly).clusters;
    EXPECT_EQ(ni_ni_fe_triplets.size(), 2);

    // Other combinations like {"Fe", "Ni", "Ni"} are impossible as no Fe has two Ni neighbors.
    AtomicCluster3Expansion fe_ni_ni(Species3AtomicSorted("Fe|Ni,Ni"));
    auto fe_ni_ni_triplets = fe_ni_ni.findAll(nl, CalculationType::ValueOnly).clusters;
    EXPECT_EQ(fe_ni_ni_triplets.size(), 0);
}

TEST(TestClusterExpansions, Cluster3Expansion) {
    // Equilateral triangle of 3 atoms, side 1.
    // Cutoff is 1.5, so all are mutually within cutoff (pure 3-body).
    Atoms atoms({{0, 0, 0}, {1, 0, 0}, {0.5, 0.866025, 0}}, {Species("Fe"), Species("Fe"), Species("Fe")});
    auto nl = NeighbourLists(atoms, 1.5);

    // There's exactly 1 unordered triplet of Fe-Fe-Fe.
    Cluster3Expansion fe_fe_fe(Species3Sorted("Fe,Fe,Fe"));
    auto fe_fe_fe_triplets = fe_fe_fe.find(nl, CalculationType::ValueOnly).clusters;
    EXPECT_EQ(fe_fe_fe_triplets.size(), 1);
}

TEST(TestClusterExpansions, ClusterPermutationModeReducedVsPermuteSameSpecies) {
    Atoms atoms({{0, 0, 0}, {1, 0, 0}, {0.5, 0.866025, 0}}, {Species("Fe"), Species("Fe"), Species("Fe")});
    auto nl = NeighbourLists(atoms, 1.5);

    AtomicCluster3Expansion reduced(Species3AtomicSorted("Fe|Fe,Fe"), ClusterPermutationMode::Reduced);
    EXPECT_EQ(reduced.getPermutationReductionFactor(), 2.0);
    auto reduced_clusters = reduced.find(0, nl, CalculationType::WithGradients);
    EXPECT_EQ(reduced_clusters.clusters.size(), 1);

    AtomicCluster3Expansion permute(Species3AtomicSorted("Fe|Fe,Fe"), ClusterPermutationMode::PermuteSameSpecies);
    EXPECT_EQ(permute.getPermutationReductionFactor(), 1.0);
    auto permute_clusters = permute.find(0, nl, CalculationType::WithGradients);
    EXPECT_EQ(permute_clusters.clusters.size(), 2);

    // The two permuted clusters share the same central atom, but swap neighbor 1 and neighbor 2.
    const auto& c0 = permute_clusters.clusters[0];
    const auto& c1 = permute_clusters.clusters[1];

    EXPECT_EQ(c0.atom_indexes[0], c1.atom_indexes[0]);
    EXPECT_EQ(c0.atom_indexes[1], c1.atom_indexes[2]);
    EXPECT_EQ(c0.atom_indexes[2], c1.atom_indexes[1]);

    EXPECT_DOUBLE_EQ(c0.separation_magnitudes[0], c1.separation_magnitudes[1]);
    EXPECT_DOUBLE_EQ(c0.separation_magnitudes[1], c1.separation_magnitudes[0]);
    EXPECT_DOUBLE_EQ(c0.separation_magnitudes[2], c1.separation_magnitudes[2]);

    // Derivatives for r12 (index 2) point in opposite directions in the permuted cluster
    const auto& d0 = permute_clusters.derivatives->at(0);
    const auto& d1 = permute_clusters.derivatives->at(1);
    EXPECT_DOUBLE_EQ((d0.val[2].direction + d1.val[2].direction).norm(), 0.0);

    // Reduction factor is always 2.0 when nodes are different species (e.g. Fe|Cu,Al)
    AtomicCluster3Expansion reduced_diff_species(Species3AtomicSorted("Fe|Al,Cu"), ClusterPermutationMode::Reduced);
    EXPECT_EQ(reduced_diff_species.getPermutationReductionFactor(), 2.0);
    AtomicCluster3Expansion reduced_diff_species2(Species3AtomicSorted("Fe|Al,Cu"),
                                                  ClusterPermutationMode::PermuteSameSpecies);
    EXPECT_EQ(reduced_diff_species.getPermutationReductionFactor(), 2.0);
}
