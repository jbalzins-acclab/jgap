#ifndef JGAP_ATOMICCLUSTER3EXPANSION_HPP
#define JGAP_ATOMICCLUSTER3EXPANSION_HPP

#include <vector>

#include "Cluster3ExpansionResult.hpp"
#include "core/CalculationType.hpp"
#include "core/atomic/neighbours/NeighbourLists.hpp"
#include "core/atomic/species/composition/Species3AtomicSorted.hpp"
#include "utils/Utils.hpp"

namespace jgap {

    class AtomicCluster3Expansion {
    public:
        static constexpr size_t ClusterPermutations = factorial(2u);

        AtomicCluster3Expansion(const Species3AtomicSorted& species_set) : species_set(species_set) {}

        void find_append(size_t atom_index, const NeighbourLists& neighbour_list,
                         Cluster3ExpansionResult& result) const {
            const Species species1 = species_set.nodes[0];
            const Species species2 = species_set.nodes[1];

            const auto sep_list1_it = neighbour_list.neighbours_per_atom[atom_index].find(species1);
            if (sep_list1_it != neighbour_list.neighbours_per_atom[atom_index].end()) {
                const auto sep_list2_it = neighbour_list.neighbours_per_atom[atom_index].find(species2);

                if (sep_list2_it != neighbour_list.neighbours_per_atom[atom_index].end()) {
                    const auto& sep_list1 = sep_list1_it->second;
                    const auto& sep_list2 = sep_list2_it->second;

                    if (species1 == species2) {
                        result.reserve(result.clusters.size() + sep_list1.size() * (sep_list2.size() - 1) / 2);
                        for (auto it1 = sep_list1.begin(); it1 != sep_list1.end(); ++it1) {
                            for (auto it2 = std::next(it1); it2 != sep_list2.end(); ++it2) {
                                std::array neigh_array{*it1, *it2};
                                result.add(atom_index, neigh_array);
                            }
                        }
                    } else {
                        result.reserve(result.clusters.size() + sep_list1.size() * sep_list2.size());
                        for (const auto& it1: sep_list1) {
                            for (const auto& it2: sep_list2) {
                                std::array neigh_array{it1, it2};
                                result.add(atom_index, neigh_array);
                            }
                        }
                    }
                }
            }
        }

        Cluster3ExpansionResult find(size_t atom_index, const NeighbourLists& neighbour_list,
                                     CalculationType calc_type) const {
            Cluster3ExpansionResult result(calc_type);
            find_append(atom_index, neighbour_list, result);
            return result;
        }

        Cluster3ExpansionResult findAll(const NeighbourLists& neighbour_list, CalculationType calc_type) const {
            Cluster3ExpansionResult result(calc_type);
            const auto it = neighbour_list.atoms_by_species.find(species_set.root);
            if (it != neighbour_list.atoms_by_species.end()) {
                for (size_t atom_index: it->second) {
                    find_append(atom_index, neighbour_list, result);
                }
            }
            return result;
        }

    private:
        Species3AtomicSorted species_set;
    };
}

#endif
