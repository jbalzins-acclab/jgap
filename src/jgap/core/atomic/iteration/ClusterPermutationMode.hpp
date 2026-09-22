#ifndef JGAP_CLUSTERPERMUTATIONMODE_HPP
#define JGAP_CLUSTERPERMUTATIONMODE_HPP

namespace jgap {
    /// @brief Indicates whether i<j or i!=j iteration should be used in cluster iteration when
    /// atomic species match at two or more indices.
    ///
    /// The main JGAP convention for cluster expansion is to provide a single cluster once,
    /// providing atoms so that node "indexes" would go in ascending order.
    /// In addition to the index of an atom in the lattice basis,
    /// the complete index implicitly includes species and
    ///
    enum class ClusterPermutationMode {
        ///
        NoNodePermutation,
        ///
        PermuteSameSpeciesNodes
    };
}

#endif
