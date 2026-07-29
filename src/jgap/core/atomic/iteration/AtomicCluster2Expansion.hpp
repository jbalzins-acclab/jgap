#ifndef JGAP_ATOMICCLUSTER2EXPANSION_HPP
#define JGAP_ATOMICCLUSTER2EXPANSION_HPP

#include <vector>

#include "jgap/core/atomic/geometry/Cluster2.hpp"
#include "jgap/core/atomic/neighbours/NeighbourLists.hpp"
#include "jgap/core/atomic/species/composition/Species2Atomic.hpp"
#include "jgap/utils/Utils.hpp"

namespace jgap {

    class AtomicCluster2Expansion {
    public:
        static constexpr Real ClusterPermutationsAvailable = 1.0;

        AtomicCluster2Expansion(const Species2Atomic& species_set) : species_set(species_set) {}

        bool forEach(size_t atom_index, const NeighbourLists& neighbour_list, Cluster2Callback auto&& callback) const {
            const auto& node_species = species_set.node;
            auto atom_neighbours = neighbour_list.neighbours_per_atom[atom_index].find(node_species);

            if (atom_neighbours == neighbour_list.neighbours_per_atom[atom_index].end()) return false;
            if (atom_neighbours->second.empty()) return false;

            for (const auto& neigh_data: atom_neighbours->second) {
                callback(
                    Cluster2{
                        .idx0 = atom_index, .idx1 = neigh_data.neighbour_index, .separation01 = neigh_data.separation
                    }
                );
            }
            return true;
        }

        bool forEach(const NeighbourLists& neighbour_list, Cluster2Callback auto&& callback) const {
            auto it = neighbour_list.atoms_by_species.find(species_set.root);
            if (it == neighbour_list.atoms_by_species.end()) return false;

            bool found = false;
            for (size_t atom_index: it->second) {
                found |= forEach(atom_index, neighbour_list, callback);
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

        Species2Atomic getSpecies() const { return species_set; }

    private:
        Species2Atomic species_set;
    };
}

#endif
