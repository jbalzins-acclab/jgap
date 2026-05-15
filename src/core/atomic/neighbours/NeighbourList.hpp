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
#include "core/atomic/geometry/Separations.hpp"

namespace jgap {

    using AtomNeighbourLists = std::map<Species, std::vector<NeighbourData>>;

    struct NeighbourList {

        static std::vector<NeighbourList> form(const std::vector<Atoms>& boxes, Real cutoff);

        template<typename... Transformers>
        static std::vector<NeighbourList> generateFor(const std::vector<Atoms>& atoms,
                                                      const Transformers&... transformers) {
            Real max_cutoff = 0.0;
            ([&]{
                max_cutoff = std::max(max_cutoff, transformers->getCutoff());
            }(), ...);
            return form(atoms, max_cutoff);
        }

        Real cutoff{};
        std::vector<AtomNeighbourLists> neighbours_per_atom;
        std::map<Species, std::vector<size_t>> atoms_by_species;

        NeighbourList() = default;
        NeighbourList(const Atoms& box, Real cutoff);

        template<size_t AtomsConnected>
        requires(AtomsConnected > 1)
        std::vector<Separations<AtomsConnected>> findAllSeparations(const SpeciesSet& species_set,
            std::optional<Real> max_distance = std::nullopt) const;

        template<size_t AtomsConnected>
        requires(AtomsConnected > 0)
        std::vector<SpeciesSet> getSpeciesSets() const;

        size_t nAtoms() const {
            return neighbours_per_atom.size();
        }

        Species speciesOf(size_t atom_index) const;

    private:
        static std::array<int, 3> findMaxRep(const Atoms& structure, Real cutoff);
    };
}

#endif
