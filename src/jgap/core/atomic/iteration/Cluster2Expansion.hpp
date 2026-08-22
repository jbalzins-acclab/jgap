#ifndef JGAP_CLUSTER2EXPANSION_HPP
#define JGAP_CLUSTER2EXPANSION_HPP

#include <optional>
#include <vector>

#include "jgap/core/atomic/geometry/Cluster2.hpp"
#include "jgap/core/atomic/neighbours/NeighbourLists.hpp"
#include "jgap/core/atomic/species/composition/Species2Atomic.hpp"
#include "jgap/core/atomic/species/composition/Species2Sorted.hpp"
#include "jgap/utils/Utils.hpp"

namespace jgap {

    class Cluster2Expansion {
    public:
        static constexpr Real ClusterPermutationsAvailable = 1.0;

        Cluster2Expansion(const Species2Atomic& species_set) : species_set1(species_set), species_set2(std::nullopt) {}

        Cluster2Expansion(const Species2Sorted& species_sorted) :
            species_set1(species_sorted.nodes[0], species_sorted.nodes[1]),
            species_set2(
                species_sorted.nodes[0] == species_sorted.nodes[1]
                    ? std::nullopt
                    : std::optional<Species2Atomic>{Species2Atomic{species_sorted.nodes[1], species_sorted.nodes[0]}}
            ) {}

        bool forEach(size_t atom_index, const NeighbourLists& neighbour_list, Cluster2Callback auto&& callback) const {
            bool found = forEach(species_set1, atom_index, neighbour_list, callback);
            if (species_set2.has_value()) {
                found |= forEach(*species_set2, atom_index, neighbour_list, callback);
            }
            return found;
        }

        bool forEach(const NeighbourLists& neighbour_list, Cluster2Callback auto&& callback) const {
            bool found = forEach(species_set1, neighbour_list, callback);
            if (species_set2.has_value()) {
                found |= forEach(*species_set2, neighbour_list, callback);
            }
            return found;
        }

        std::vector<Cluster2> expand(size_t atom_index, const NeighbourLists& neighbour_list) const {
            std::vector<Cluster2> result;
            forEach(atom_index, neighbour_list, [&](const Cluster2& cluster) { result.emplace_back(cluster); });
            return result;
        }

        std::vector<Cluster2> expand(const NeighbourLists& neighbour_list) const {
            std::vector<Cluster2> result;
            forEach(neighbour_list, [&](const Cluster2& cluster) { result.emplace_back(cluster); });
            return result;
        }

        Species2Atomic getSpecies() const { return species_set1; }
        const std::optional<Species2Atomic>& getSecondarySpecies() const { return species_set2; }

    private:
        static bool forEach(
            const Species2Atomic& sp, size_t atom_index, const NeighbourLists& neighbour_list,
            Cluster2Callback auto&& callback
        ) {
            const auto& node_species = sp.node;
            auto atom_neighbours = neighbour_list.neighbours_per_atom[atom_index].find(node_species);

            if (atom_neighbours == neighbour_list.neighbours_per_atom[atom_index].end()) return false;
            if (atom_neighbours->second.empty()) return false;

            for (const auto& neigh_data: atom_neighbours->second) {
                callback(
                    Cluster2{
                        .idx0 = atom_index,
                        .idx1 = neigh_data.neighbour_index,
                        .separation01 = neigh_data.separation,
                    }
                );
            }
            return true;
        }

        static bool forEach(
            const Species2Atomic& sp, const NeighbourLists& neighbour_list, Cluster2Callback auto&& callback
        ) {
            auto it = neighbour_list.atoms_by_species.find(sp.root);
            if (it == neighbour_list.atoms_by_species.end()) return false;

            bool found = false;
            for (size_t atom_index: it->second) {
                found |= forEach(sp, atom_index, neighbour_list, callback);
            }
            return found;
        }

        Species2Atomic species_set1;
        std::optional<Species2Atomic> species_set2;
    };
}

#endif
