#ifndef JGAP_NEIGHBOURLIST_HPP
#define JGAP_NEIGHBOURLIST_HPP

#include <map>
#include <vector>
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

    using AtomNeighbourLists = std::map<Species, std::vector<NeighbourData>>;

    struct NeighbourList {

        static std::vector<NeighbourList> form(const std::vector<Atoms>& boxes, Real cutoff);

        template<typename... Transformers>
        static std::vector<NeighbourList> generateFor(const std::vector<Atoms>& atoms,
                                                      const Transformers&... transformers) {
            Cutoffs combined_cutoffs;
            (..., (combined_cutoffs += transformers.getCutoffs()));
            return form(atoms, combined_cutoffs.maxOverall());
        }

        Real cutoff{};
        std::vector<AtomNeighbourLists> neighbours_per_atom;
        std::map<Species, std::vector<size_t>> atoms_by_species;

        NeighbourList() = default;
        NeighbourList(const Atoms& box, Real cutoff);

        template<size_t N, ClusterTypes ClusterType>
        requires(N > 1 && N <= 3)
        std::vector<Cluster<N>> findAllClusters(const SpeciesSet<N, ClusterType>& species_set)
            const;

        template<size_t N, ClusterTypes ClusterType>
        requires(N > 1 && N <= 3)
        std::vector<SpeciesSet<N, ClusterType>> getSpeciesSets() const;

        template<size_t N>
        requires(N > 1 && N <= 3)
        std::vector<SpeciesSet<N, HasCentralAtom>> getSpeciesSets(Species central_atom_species) const;

        size_t nAtoms() const {
            return neighbours_per_atom.size();
        }

        Species speciesOf(size_t atom_index) const;

    private:
        static std::array<int, 3> findMaxRep(const Atoms& structure, Real cutoff);
    };
}

#endif