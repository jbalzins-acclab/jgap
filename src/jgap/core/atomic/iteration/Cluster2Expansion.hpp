#ifndef JGAP_CLUSTER2EXPANSION_HPP
#define JGAP_CLUSTER2EXPANSION_HPP

#include <vector>

#include "jgap/core/atomic/geometry/Cluster2.hpp"
#include "jgap/core/atomic/neighbours/NeighbourLists.hpp"
#include "jgap/core/atomic/species/composition/Species2Sorted.hpp"
#include "jgap/utils/Utils.hpp"

namespace jgap {

    class Cluster2Expansion {
    public:
        static constexpr size_t ClusterPermutationsAvailable = factorial(2);

        Cluster2Expansion(const Species2Sorted& species_set) : species_set(species_set) {}

        bool forEach(const NeighbourLists& neighbour_list, Cluster2Callback auto&& callback) const {
            for (auto& species: species_set.nodes) {
                if (!neighbour_list.atoms_by_species.contains(species)) return false;
            }

            const auto& species1 = species_set.nodes[0];
            const auto& species2 = species_set.nodes[1];
            const auto& atoms = neighbour_list.atoms_by_species.at(species1);
            if (atoms.empty()) return false;

            for (size_t i: atoms) {
                auto atom_neighbours = neighbour_list.neighbours_per_atom[i].find(species2);

                if (atom_neighbours == neighbour_list.neighbours_per_atom[i].end()) continue;

                for (auto neigh_data: atom_neighbours->second) {
                    if (species1 == species2) {
                        if (neigh_data.neighbour_index == i) [[unlikely]] {
                            if (neigh_data.separation.direction.isPositive()) {
                                callback(
                                    Cluster2{
                                        .idx0 = i,
                                        .idx1 = neigh_data.neighbour_index,
                                        .separation01 = neigh_data.separation
                                    }
                                );
                            }
                        } else if (neigh_data.neighbour_index > i) {
                            callback(
                                Cluster2{
                                    .idx0 = i, .idx1 = neigh_data.neighbour_index, .separation01 = neigh_data.separation
                                }
                            );
                        }
                    } else {
                        callback(
                            Cluster2{
                                .idx0 = i, .idx1 = neigh_data.neighbour_index, .separation01 = neigh_data.separation
                            }
                        );
                    }
                }
            }
            return true;
        }

        std::vector<Cluster2> expand(const NeighbourLists& neighbour_list) const {
            std::vector<Cluster2> result;

            forEach(neighbour_list, [&](const Cluster2& cluster) { result.emplace_back(cluster); });

            return result;
        }

        Species2Sorted getSpecies() const { return species_set; }

    private:
        Species2Sorted species_set;
    };
}

#endif
