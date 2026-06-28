#ifndef JGAP_CLUSTERFINDER_HPP
#define JGAP_CLUSTERFINDER_HPP
#include <vector>

#include "core/Real.hpp"
#include "core/atomic/geometry/Cluster.hpp"
#include "core/atomic/geometry/ClusterSymmetry.hpp"
#include "core/atomic/neighbours/NeighbourList.hpp"
#include "core/atomic/species/SpeciesSet.hpp"

namespace jgap {

    template<size_t ClusterSize, ClusterSymmetry ClusterSym>
    class ClusterFinder {
    public:
        static Real getRepetitionReductionFactor() {
            if constexpr (ClusterSym == IndexSensitive) return 0;
            return factorial(ClusterSize - (ClusterSym == NodeSymmetric));
        }

        /*
        template<CalculationType CalcType>
        static std::vector<Cluster<2, CalcType>> find(const NeighbourList& neighbour_list,
                                                      const SpeciesSet<2, FullSymmetry>& species_set) {

            for (auto& species: species_set.getNodes()) {
                if (!neighbour_list.atoms_by_species.contains(species)) return {};
            }

            std::vector<Cluster<2, CalcType>> result;

            const auto& species1 = species_set.getNodes()[0];
            const auto& species2 = species_set.getNodes()[1];

            for (size_t i: neighbour_list.atoms_by_species.at(species1)) {
                auto atom_neighbours = neighbour_list.neighbours_per_atom[i].find(species2);

                if (atom_neighbours == neighbour_list.neighbours_per_atom[i].end()) continue;

                for (auto neigh_data: atom_neighbours->second) {

                    if (neigh_data.neighbour_index == i) [[unlikely]] {
                        if (neigh_data.separation.derivatives.direction.x < 0.0) continue;
                        if (neigh_data.separation.derivatives.direction.x > 0.0) {
                            result.emplace_back(i, std::array{neigh_data});
                            continue;
                        }
                        if (neigh_data.separation.derivatives.direction.y < 0.0) continue;
                        if (neigh_data.separation.derivatives.direction.y > 0.0) {
                            result.emplace_back(i, std::array{neigh_data});
                            continue;
                        }
                        if (neigh_data.separation.derivatives.direction.z < 0.0) continue;
                        result.emplace_back(i, std::array{neigh_data});
                    } else if (neigh_data.neighbour_index > i) {
                        result.emplace_back(i, std::array{neigh_data});
                    }
                }
            }

            return result;
        }

        template<CalculationType CalcType, ClusterSymmetry ClusterSym>
        requires(ClusterSym == NodeSymmetric || ClusterSym == IndexSensitive)
        static std::vector<Cluster<2, CalcType>> find(const NeighbourList& neighbour_list,
                                                      const SpeciesSet<2, ClusterSym>& species_set) {

            if (!neighbour_list.atoms_by_species.contains(species_set.getRoot())) return {};

            for (auto& species: species_set.getNodes()) {
                if (!neighbour_list.atoms_by_species.contains(species)) return {};
            }

            std::vector<Cluster<2, CalcType>> result;

            const auto& central_species = species_set.getRoot();
            const auto& node_species = species_set.getNodes()[0];

            for (size_t i: neighbour_list.atoms_by_species.at(central_species)) {
                auto atom_neighbours = neighbour_list.neighbours_per_atom[i].find(node_species);

                if (atom_neighbours == neighbour_list.neighbours_per_atom[i].end()) continue;

                for (auto neigh_data: atom_neighbours->second) {
                    result.emplace_back(i, std::array{neigh_data});
                }
            }

            return result;
        }

        template<CalculationType CalcType>
        static std::vector<Cluster<3, CalcType>> find(const NeighbourList& neighbour_list,
                                                      const SpeciesSet<3, NodeSymmetric>& species_set) {


            if (!neighbour_list.atoms_by_species.contains(species_set.getRoot())) return {};

            for (auto& species: species_set.getNodes()) {
                if (!neighbour_list.atoms_by_species.contains(species)) return {};
            }

            std::vector<Cluster<3, CalcType>> result;

            const auto& central_species = species_set.getRoot();
            const auto& species1 = species_set.getNodes()[0];
            const auto& species2 = species_set.getNodes()[1];

            for (const size_t i: neighbour_list.atoms_by_species.at(central_species)) {

                const auto sep_list1_it = neighbour_list.neighbours_per_atom[i].find(species1);
                if (sep_list1_it == neighbour_list.neighbours_per_atom[i].end()) continue;

                const auto sep_list2_it = neighbour_list.neighbours_per_atom[i].find(species2);
                if (sep_list2_it == neighbour_list.neighbours_per_atom[i].end()) continue;

                const auto& sep_list1 = sep_list1_it->second;
                const auto& sep_list2 = sep_list2_it->second;

                if (species1 == species2) {
                    for (auto it1 = sep_list1.begin(); it1 != sep_list1.end(); ++it1) {
                        for (auto it2 = std::next(it1); it2 != sep_list2.end(); ++it2) {
                            result.emplace_back(i, std::array{*it1, *it2});
                        }
                    }
                } else {
                    for (const auto& it1 : sep_list1) {
                        for (const auto& it2 : sep_list2) {
                            result.emplace_back(i, std::array{it1, it2});
                        }
                    }
                }
            }

            return result;
        }

        template<CalculationType CalcType>
        static std::vector<Cluster<3, CalcType>> find(const NeighbourList& neighbour_list,
                                                      const SpeciesSet<3, IndexSensitive>& species_set) {


            if (!neighbour_list.atoms_by_species.contains(species_set.getRoot())) return {};

            for (auto& species: species_set.getNodes()) {
                if (!neighbour_list.atoms_by_species.contains(species)) return {};
            }

            std::vector<Cluster<3, CalcType>> result;

            const auto& central_species = species_set.getRoot();
            const auto& species1 = species_set.getNodes()[0];
            const auto& species2 = species_set.getNodes()[1];

            for (const size_t i: neighbour_list.atoms_by_species.at(central_species)) {
                if (species1 == species2) {

                    const auto sep_list_it = neighbour_list.neighbours_per_atom[i].find(species1);
                    if (sep_list_it == neighbour_list.neighbours_per_atom[i].end()) [[unlikely]] continue;

                    const auto& sep_list1 = sep_list_it->second;

                    for (auto it1 = sep_list1.begin(); it1 != sep_list1.end(); ++it1) {
                        for (auto it2 = sep_list1.begin(); it2 != sep_list1.end(); ++it2) {
                            if (it1 == it2) [[unlikely]] continue;
                            result.emplace_back(i, std::array{*it1, *it2});
                        }
                    }

                } else {

                    const auto sep_list1_it = neighbour_list.neighbours_per_atom[i].find(species1);
                    if (sep_list1_it == neighbour_list.neighbours_per_atom[i].end()) [[unlikely]] continue;

                    const auto sep_list2_it = neighbour_list.neighbours_per_atom[i].find(species2);
                    if (sep_list2_it == neighbour_list.neighbours_per_atom[i].end()) [[unlikely]] continue;

                    const auto& sep_list1 = sep_list1_it->second;
                    const auto& sep_list2 = sep_list2_it->second;

                    for (const auto& it1: sep_list1) {
                        for (const auto& it2: sep_list2) {
                            result.emplace_back(i, std::array{it1, it2});
                        }
                    }
                }
            }

            return result;
        }*/
    };
}

#endif
