#ifndef JGAP_CLUSTER2EXPANSION_HPP
#define JGAP_CLUSTER2EXPANSION_HPP

#include <vector>

#include "Cluster2ExpansionResult.hpp"
#include "core/CalculationType.hpp"
#include "core/atomic/neighbours/NeighbourLists.hpp"
#include "core/atomic/species/composition/Species2Sorted.hpp"
#include "utils/Utils.hpp"

namespace jgap {

    class Cluster2Expansion {
    public:
        static constexpr size_t ClusterPermutations = factorial(2);

        Cluster2Expansion(const Species2Sorted &species_set) : species_set(species_set) {}

        Cluster2ExpansionResult find(const NeighbourLists &neighbour_list, CalculationType calc_type) const {
            Cluster2ExpansionResult result(calc_type);

            for (auto &species: species_set.nodes) {
                if (!neighbour_list.atoms_by_species.contains(species)) return result;
            }

            const auto &species1 = species_set.nodes[0];
            const auto &species2 = species_set.nodes[1];

            for (size_t i: neighbour_list.atoms_by_species.at(species1)) {
                auto atom_neighbours = neighbour_list.neighbours_per_atom[i].find(species2);

                if (atom_neighbours == neighbour_list.neighbours_per_atom[i].end()) continue;

                result.reserve(result.clusters.size() + atom_neighbours->second.size());

                for (auto neigh_data: atom_neighbours->second) {
                    if (neigh_data.neighbour_index == i) [[unlikely]] {
                        if (neigh_data.separation.derivatives.direction.isPositive()) {
                            std::array neigh_array{neigh_data};
                            result.add(i, neigh_array);
                        }
                    } else if (neigh_data.neighbour_index > i) {
                        std::array neigh_array{neigh_data};
                        result.add(i, neigh_array);
                    }
                }
            }

            return result;
        }

        const Species2Sorted& getSpecies() const { return species_set; }

    private:
        Species2Sorted species_set;
    };
}

#endif
