#ifndef JGAP_NEIGHBOURLIST_HPP
#define JGAP_NEIGHBOURLIST_HPP

#include <array>
#include <map>
#include <vector>

#include "NeighbourData.hpp"
#include "jgap/core/atomic/Atoms.hpp"
#include "jgap/core/atomic/species/Species.hpp"

namespace jgap {

    /// \brief Formation, storage, and manipulation of per-atom neighbour lists.
    ///
    /// Stores per-atom per-species neighbour lists. <br>
    ///
    /// \warning Formed during construction. <br>
    ///
    /// \note Grouping neighbour lists per species is crucial for simplifying species-dependent routing logic
    /// requiring iteration.
    ///
    struct NeighbourLists {
        using AtomNeighbourLists = std::map<Species, std::vector<NeighbourData>>;

        static std::vector<NeighbourLists> form(const std::vector<Atoms> &boxes, Real cutoff);

        std::vector<AtomNeighbourLists> neighbours_per_atom;
        std::map<Species, std::vector<size_t>> atoms_by_species;

        NeighbourLists(const Atoms &box, Real cutoff);

        size_t nAtoms() const { return neighbours_per_atom.size(); }
        Real getCutoff() const { return cutoff; }

    private:
        static std::array<int, 3> findMaxRep(const Atoms &structure, Real cutoff);

        Real cutoff{};
    };
}

#endif
