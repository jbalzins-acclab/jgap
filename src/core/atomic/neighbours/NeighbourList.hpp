#ifndef JGAP_NEIGHBOURLIST_HPP
#define JGAP_NEIGHBOURLIST_HPP

#include <map>
#include <vector>
#include <functional>
#include <tuple>
#include <optional>
#include <array>

#include "NeighbourData.hpp"
#include "core/atomic/species/Species.hpp"
#include "core/atomic/species/SpeciesSet.hpp"
#include "core/atomic/Box.hpp"
#include "core/atomic/geometry/Separations.hpp"

namespace jgap {

    using AtomNeighbourLists = std::map<Species, std::vector<NeighbourData>>;

    struct NeighbourList {

        static std::vector<NeighbourList> findAllNeighbours(const std::vector<Box>& boxes, double cutoff);

        double cutoff;
        std::vector<AtomNeighbourLists> neighbours_per_atom;
        std::map<Species, std::vector<size_t>> atoms_by_species;

        NeighbourList() = default;
        NeighbourList(const Box& box, double cutoff);

        template<size_t N>
        requires(N > 1)
        std::vector<Separations<N>> findAllSeparations(const SpeciesSet& species_set,
            std::optional<double> max_distance = std::nullopt) const;

        template<size_t N>
        requires(N > 0)
        std::vector<SpeciesSet> getSpeciesSets() const;

        size_t nAtoms() const {
            return neighbours_per_atom.size();
        }

    private:
        static std::array<int, 3> findMaxRep(const Box& structure, double cutoff);
    };
}

#endif
