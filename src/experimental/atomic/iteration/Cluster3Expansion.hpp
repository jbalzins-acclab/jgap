#ifndef JGAP_CLUSTER3EXPANSION_HPP
#define JGAP_CLUSTER3EXPANSION_HPP

#include <vector>

#include "core/CalculationType.hpp"
#include "core/atomic/iteration/Cluster3ExpansionResult.hpp"
#include "core/atomic/neighbours/NeighbourLists.hpp"
#include "experimental/atomic/species/composition/Species3Sorted.hpp"
#include "utils/Utils.hpp"

namespace jgap {

    class Cluster3Expansion {
    public:
        static constexpr size_t ClusterPermutations = factorial(3u);

        Cluster3Expansion(const Species3Sorted& species_set) : species_set(species_set) {}

        Cluster3ExpansionResult find(const NeighbourLists& neighbour_list, CalculationType calc_type) const {
            Cluster3ExpansionResult result(calc_type);

            for (auto& species: species_set.nodes) {
                if (!neighbour_list.atoms_by_species.contains(species)) return result;
            }

            const auto& species1 = species_set.nodes[0];
            const auto& species2 = species_set.nodes[1];
            const auto& species3 = species_set.nodes[2];

            Real cutoff = neighbour_list.getCutoff();

            for (size_t i: neighbour_list.atoms_by_species.at(species1)) {
                auto atom_neighbours_j = neighbour_list.neighbours_per_atom[i].find(species2);
                if (atom_neighbours_j == neighbour_list.neighbours_per_atom[i].end()) continue;

                auto atom_neighbours_k = neighbour_list.neighbours_per_atom[i].find(species3);
                if (atom_neighbours_k == neighbour_list.neighbours_per_atom[i].end()) continue;

                const auto& list_j = atom_neighbours_j->second;
                const auto& list_k = atom_neighbours_k->second;

                for (const auto& neigh_j: list_j) {
                    size_t j = neigh_j.neighbour_index;

                    if (species1 == species2) {
                        if (j < i) continue;
                        if (j == i && !neigh_j.separation.vec().isPositive()) continue;
                    }

                    for (const auto& neigh_k: list_k) {
                        size_t k = neigh_k.neighbour_index;

                        if (species2 == species3) {
                            if (k < j) continue;
                            if (k == j) {
                                Vector3 r_jk = neigh_k.separation.vec() - neigh_j.separation.vec();
                                if (!r_jk.isPositive()) continue;
                            }
                        }

                        // Also check that j-k separation is < cutoff
                        Vector3 r_jk = neigh_k.separation.vec() - neigh_j.separation.vec();
                        if (r_jk.norm() >= cutoff) continue;

                        std::array neigh_array{neigh_j, neigh_k};
                        result.add(i, neigh_array);
                    }
                }
            }

            return result;
        }

        const Species3Sorted& getSpecies() const { return species_set; }

    private:
        Species3Sorted species_set;
    };
}

#endif
