#ifndef JGAP_CLUSTER3EXPANSION_HPP
#define JGAP_CLUSTER3EXPANSION_HPP

#include <vector>

#include "jgap/core/atomic/geometry/Cluster3.hpp"
#include "jgap/core/atomic/neighbours/NeighbourLists.hpp"
#include "jgap/experimental/atomic/species/composition/Species3Sorted.hpp"
#include "jgap/utils/Utils.hpp"

namespace jgap {

    class Cluster3Expansion {
    public:
        static constexpr Real ClusterPermutationsAvailable = factorial(3u);

        Cluster3Expansion(const Species3Sorted& species_set) : species_set(species_set) {}

        bool forEach(const NeighbourLists& neighbour_list, Cluster3Callback auto&& callback) const {
            for (auto& species: species_set.nodes) {
                if (!neighbour_list.atoms_by_species.contains(species)) return false;
            }

            const auto& species1 = species_set.nodes[0];
            const auto& species2 = species_set.nodes[1];
            const auto& species3 = species_set.nodes[2];

            Real cutoff = neighbour_list.getCutoff();
            bool found = false;

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
                        if (j == i && !neigh_j.separation.direction.isPositive()) continue;
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

                        Separation sep12(neigh_j.separation.vec(), neigh_k.separation.vec());
                        if (sep12.magnitude >= cutoff) continue;

                        found = true;
                        callback(
                            Cluster3{
                                .atom_indexes = {i, j, k},
                                .separations = {neigh_j.separation, neigh_k.separation, sep12}
                            }
                        );
                    }
                }
            }
            return found;
        }

        std::vector<Cluster3> expand(const NeighbourLists& neighbour_list) const {
            std::vector<Cluster3> result;

            forEach(neighbour_list, [&](const Cluster3& cluster) { result.emplace_back(cluster); });

            return result;
        }

        Species3Sorted getSpecies() const { return species_set; }

    private:
        Species3Sorted species_set;
    };
}

#endif
