#ifndef JGAP_CLUSTER3ATOMICEXPANSION_HPP
#define JGAP_CLUSTER3ATOMICEXPANSION_HPP

#include <vector>

#include "ClusterPermutationMode.hpp"
#include "jgap/core/atomic/geometry/Cluster3.hpp"
#include "jgap/core/atomic/neighbours/NeighbourLists.hpp"
#include "jgap/core/atomic/species/composition/Species3AtomicSorted.hpp"
#include "jgap/utils/Utils.hpp"

namespace jgap {

    class Cluster3AtomicExpansion {
    public:
        static constexpr Real ClusterPermutationsAvailable = factorial(2u);

        Cluster3AtomicExpansion(
            const Species3AtomicSorted& species_set, ClusterPermutationMode mode = ClusterPermutationMode::Reduced
        ) :
            species_set(species_set),
            mode(mode),
            factor(
                (species_set.nodes[0] != species_set.nodes[1] || mode == ClusterPermutationMode::Reduced)
                    ? static_cast<Real>(ClusterPermutationsAvailable)
                    : 1.0
            ) {}

        Real getPermutationReductionFactor() const { return factor; }

        ClusterPermutationMode getMode() const { return mode; }

        bool forEach(size_t atom_index, const NeighbourLists& neighbour_list, Cluster3Callback auto&& callback) const {
            const Species species1 = species_set.nodes[0];
            const Species species2 = species_set.nodes[1];

            const auto sep_list1_it = neighbour_list.neighbours_per_atom[atom_index].find(species1);
            if (sep_list1_it == neighbour_list.neighbours_per_atom[atom_index].end()) return false;

            const auto& sep_list1 = sep_list1_it->second;
            bool found;

            if (species1 == species2) {
                found = sep_list1.size() >= 2;

                for (auto it1 = sep_list1.begin(); it1 != sep_list1.end(); it1++) {
                    for (auto it2 = std::next(it1); it2 != sep_list1.end(); it2++) {
                        found = true;
                        callback(
                            Cluster3{
                                .atom_indexes = {atom_index, it1->neighbour_index, it2->neighbour_index},
                                .separations = {
                                    it1->separation, it2->separation,
                                    Separation(it1->separation.vec(), it2->separation.vec())
                                }
                            }
                        );
                        if (mode != ClusterPermutationMode::Reduced) {
                            callback(
                                Cluster3{
                                    .atom_indexes = {atom_index, it2->neighbour_index, it1->neighbour_index},
                                    .separations = {
                                        it2->separation, it1->separation,
                                        Separation(it2->separation.vec(), it1->separation.vec())
                                    }
                                }
                            );
                        }
                    }
                }
            } else {

                const auto sep_list2_it = neighbour_list.neighbours_per_atom[atom_index].find(species2);
                if (sep_list2_it == neighbour_list.neighbours_per_atom[atom_index].end()) return false;

                const auto& sep_list2 = sep_list2_it->second;

                found = !sep_list1.empty() && !sep_list2.empty();

                for (const auto& it1: sep_list1) {
                    for (const auto& it2: sep_list2) {
                        found = true;
                        callback(
                            Cluster3{
                                .atom_indexes = {atom_index, it1.neighbour_index, it2.neighbour_index},
                                .separations = {
                                    it1.separation, it2.separation,
                                    Separation(it1.separation.vec(), it2.separation.vec())
                                }
                            }
                        );
                    }
                }
            }
            return found;
        }

        bool forEach(const NeighbourLists& neighbour_list, Cluster3Callback auto&& callback) const {
            const auto it = neighbour_list.atoms_by_species.find(species_set.root);
            if (it == neighbour_list.atoms_by_species.end()) return false;

            bool found = false;
            for (size_t atom_index: it->second) {
                found |= forEach(atom_index, neighbour_list, callback);
            }
            return found;
        }

        std::vector<Cluster3> expand(size_t atom_index, const NeighbourLists& neighbour_list) const {
            std::vector<Cluster3> result;
            forEach(atom_index, neighbour_list, [&](const Cluster3& cluster) { result.emplace_back(cluster); });
            return result;
        }

        std::vector<Cluster3> expand(const NeighbourLists& neighbour_list) const {
            std::vector<Cluster3> result;
            forEach(neighbour_list, [&](const Cluster3& cluster) { result.emplace_back(cluster); });
            return result;
        }

        Species3AtomicSorted getSpecies() const { return species_set; }

    private:
        Species3AtomicSorted species_set;
        ClusterPermutationMode mode;
        Real factor;
    };
}

#endif
