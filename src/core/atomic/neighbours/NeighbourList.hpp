#ifndef JGAP_NEIGHBOURLIST_HPP
#define JGAP_NEIGHBOURLIST_HPP

#include <map>
#include <vector>
#include <set>
#include <functional>
#include <tuple>
#include <optional>
#include <array>
#include <memory>

#include "NeighbourData.hpp"
#include "core/atomic/species/Species.hpp"
#include "core/atomic/species/SpeciesSet.hpp"
#include "core/atomic/Atoms.hpp"
#include "core/atomic/geometry/Cluster.hpp"
#include "../../potentials/Cutoffs.hpp"

namespace jgap {

    /// @brief Formation, storage, and manipulation of per-atom neighbour lists.
    ///
    /// Stores per-atom per-species neighbour lists (\see NeighbourList#AtomNeighbourLists). <br>
    /// Formed during construction. <br>
    ///
    /// @note Grouping neighbour lists per species is crucial for simplifying species-dependent routing logic
    /// requiring iteration.
    ///
    struct NeighbourList {
        using AtomNeighbourLists = std::map<Species, std::vector<NeighbourData>>;

        static std::vector<NeighbourList> form(const std::vector<Atoms>& boxes, Real cutoff);

        std::vector<AtomNeighbourLists> neighbours_per_atom;
        std::map<Species, std::vector<size_t>> atoms_by_species;

        NeighbourList(const Atoms& box, Real cutoff);

        /// @brief Finds all combinations
        ///
        /// Iterating through N-body clusters in a periodic system in an efficient way poses a challenge
        /// because of the possibility of an atom's periodic image being within the relevant cutoff boundary.
        /// Re-implementing such iteration leads to mismatches between conventions,
        /// and, hence, confusion related to that.
        /// The aim of the #findAllClusters method is to standardize such iteration to avoid such confusion.
        ///
        /// Implemented for clusters of given size/symmetry as: <br>
        /// SpeciesSet<2, FullSymmetry>: |r_ij| for all i<j iteration for non-periodic systems,
        /// and it's periodic generalization; <br>
        /// SpeciesSet<2, HasCentralAtom>: |r_ij| for each neighbour j of each atom i in the system; <br>
        /// SpeciesSet<3, HasCentralAtom>: {|r_ij|, |r_ik|, |r_jk|} for each neighbour pair j<k
        /// of each atom i in the system; <br>
        ///
        ///
        /// \tparam CalcType Indicates whether derivatives should be provided in the cluster.
        /// \tparam N Size of clusters.
        /// \tparam ClusterSym Defines whether the atom at clusters index 0 should be treated
        /// as the special "root" atom which "owns" the corresponding cluster.
        /// \param species_set Defines species of which corresponding cluster atom should be.
        /// \return Clusters corresponding to all combinations(!NOT PERMUTATIONS!) of atoms
        /// of the given species set.
        ///
        /// \note Should be preferred to custom iteration to avoid any concept mismatch.
        //template<CalculationType CalcType = WithGradients, size_t N, ClusterSymmetry ClusterSym>
        //requires(N > 1 && N <= 3)
        //std::vector<Cluster<N, CalcType>> findAllClusters(const SpeciesSet<N, ClusterSym>& species_set) const;

        size_t nAtoms() const { return neighbours_per_atom.size(); }
        Real getCutoff() const { return cutoff; }

    private:
        static std::array<int, 3> findMaxRep(const Atoms& structure, Real cutoff);

        Real cutoff{};
    };
}

#endif