#include "NeighbourList.hpp"
#include <tbb/parallel_for.h>
#include <iostream>
#include <algorithm>
#include <numeric>
#include <ranges>
#include <set>
#include "io/log/CurrentLogger.hpp"

namespace jgap {

    Species NeighbourList::speciesOf(size_t atom_index) const {
        for (const auto& [species, atom_indices]: atoms_by_species) {
            for (const size_t index: atom_indices) {
                if (index == atom_index) {
                    return species;
                }
            }
        }
        throw std::runtime_error("Atom index not found in any species.");
    }

    std::array<int, 3> NeighbourList::findMaxRep(const Atoms& structure, const Real cutoff) {
        auto pbc = structure.lookupPbc();
        auto lattice_opt = structure.lookupLattice();

        if (!lattice_opt) {
            if (pbc[0] || pbc[1] || pbc[2]) {
                JGAP_LOG_AND_THROW("PBC is true but no lattice is provided.");
            }
            return {0, 0, 0};
        }
        const auto& lattice = *lattice_opt;

        const Vector3 side1 = lattice.a,
                      side2 = lattice.b,
                      side3 = lattice.c;

        std::array maxRep = {
            pbc[0] ? static_cast<int>(cutoff / side1.norm() + 1) : 0,
            pbc[1] ? static_cast<int>(cutoff / side2.norm() + 1) : 0,
            pbc[2] ? static_cast<int>(cutoff / side3.norm() + 1) : 0
        };

        // triclinic
        if (abs(side1.dot(side2)) > 1e-6 || abs(side1.dot(side3)) > 1e-6 || abs(side2.dot(side3)) > 1e-6) {
            if (pbc[0]) maxRep[0] = static_cast<int>(cutoff / side1.aproject(side2, side3)) + 1;
            if (pbc[1]) maxRep[1] = static_cast<int>(cutoff / side2.aproject(side1, side3)) + 1;
            if (pbc[2]) maxRep[2] = static_cast<int>(cutoff / side3.aproject(side2, side1)) + 1;
        }

        return maxRep;
    }

    NeighbourList::NeighbourList(const Atoms& box, Real cutoff) : cutoff(cutoff) {
        const auto max_rep = findMaxRep(box, cutoff);
        auto lattice_opt = box.lookupLattice();
        auto& species = box.lookupSpecies();
        auto& positions = box.lookupPositions();

        neighbours_per_atom.resize(box.nAtoms());

        for (size_t i = 0; i < box.nAtoms(); i++) {
            atoms_by_species[species[i]].push_back(i);
        }

        for (size_t i = 0; i < box.nAtoms(); i++) {
            auto& neighbours_i = neighbours_per_atom[i];

            for (size_t j = 0; j < box.nAtoms(); j++) {
                for (int rep0 = -max_rep[0]; rep0 <= max_rep[0]; rep0++) {
                    for (int rep1 = -max_rep[1]; rep1 <= max_rep[1]; rep1++) {
                        for (int rep2 = -max_rep[2]; rep2 <= max_rep[2]; rep2++) {

                            Vector3 offset = {0.0, 0.0, 0.0};
                            if (lattice_opt) {
                                const auto&[a, b, c] = *lattice_opt;
                                offset = a * rep0
                                       + b * rep1
                                       + c * rep2;
                            }

                            if (i == j && rep0 == 0 && rep1 == 0 && rep2 == 0) continue;

                            auto pos_j_plus_offset = positions[j] + offset;
                            auto separation_ij = Separation(positions[i], pos_j_plus_offset);

                            if (separation_ij.magnitude <= cutoff) {
                                neighbours_i[species[j]].emplace_back(j, separation_ij);
                            }
                        }
                    }
                }
            }
        }
    }

    std::vector<NeighbourList> NeighbourList::form(const std::vector<Atoms>& boxes, Real cutoff) {
        std::vector<NeighbourList> result(boxes.size());

        tbb::parallel_for(static_cast<size_t>(0), boxes.size(), [&](size_t i) {
            result[i] = NeighbourList(boxes[i], cutoff);
        });

        return result;
    }

    // CalcType first since template deduction works well for others
    template<CalculationType CalcType, size_t N, ClusterSymmetry ClusterType>
    requires(N > 1 && N <= 3)
    std::vector<Cluster<N, CalcType>> NeighbourList::findAllClusters(
        const SpeciesSet<N, ClusterType> &species_set) const {

        if constexpr (N > 3) {
            JGAP_LOG_AND_THROW("Inter-separations finding for N > 3 is not implemented");
        }

        if constexpr (ClusterType == HasCentralAtom) {
            if (!atoms_by_species.contains(species_set.getRoot())) return {};
        }
        for (auto& species: species_set.getNodes()) {
            if (!atoms_by_species.contains(species)) return {};
        }

        if constexpr (N == 2 && ClusterType == HasCentralAtom) {

            std::vector<Cluster<2, CalcType>> result;

            const auto& central_species = species_set.getRoot();
            const auto& node_species = species_set.getNodes()[0];

            for (size_t i: atoms_by_species.at(central_species)) {
                auto atom_neighbours = neighbours_per_atom[i].find(node_species);

                if (atom_neighbours == neighbours_per_atom[i].end()) continue;

                for (auto neigh_data: atom_neighbours->second) {
                    result.emplace_back(i, std::array{neigh_data});
                }
            }

            return result;

        } else if constexpr (N == 2 && ClusterType == Symmetric) {

            std::vector<Cluster<2, CalcType>> result;

            const auto& species1 = species_set.getNodes()[0];
            const auto& species2 = species_set.getNodes()[1];

            // WARN!
            // Do not try to optimize e.g. via ignoring j < i unless it's truly performance-critical
            // (which I don't expect it to be for pair potential).
            // This would require atoms that interact periodically with themselves
            // to be treated separately down the line.

            for (size_t i: atoms_by_species.at(species1)) {
                auto atom_neighbours = neighbours_per_atom[i].find(species2);

                if (atom_neighbours == neighbours_per_atom[i].end()) continue;

                for (auto neigh_data: atom_neighbours->second) {
                    result.emplace_back(i, std::array{neigh_data});
                }
            }

            if (species1 != species2) {

                // To enforce consistency with same-species iteration:
                // two clusters formed per pair.

                for (const size_t i: atoms_by_species.at(species2)) {
                    auto atom_neighbours = neighbours_per_atom[i].find(species1);

                    if (atom_neighbours == neighbours_per_atom[i].end()) continue;

                    for (auto neigh_data: atom_neighbours->second) {
                        result.emplace_back(i, std::array{neigh_data});
                    }
                }
            }

            return result;

        } else if constexpr (N == 3 && ClusterType == HasCentralAtom) {

            std::vector<Cluster<3, CalcType>> result;

            const auto& central_species = species_set.getRoot();
            const auto& species1 = species_set.getNodes()[0];
            const auto& species2 = species_set.getNodes()[1];

            for (const size_t i: atoms_by_species.at(central_species)) {

                const auto sep_list1_it = neighbours_per_atom[i].find(species1);
                if (sep_list1_it == neighbours_per_atom[i].end()) continue;

                const auto sep_list2_it = neighbours_per_atom[i].find(species2);
                if (sep_list2_it == neighbours_per_atom[i].end()) continue;

                const auto& sep_list1 = sep_list1_it->second;
                const auto& sep_list2 = sep_list2_it->second;

                for (auto it1 = sep_list1.begin(); it1 != sep_list1.end(); ++it1) {
                    auto it2_start = sep_list2.begin();
                    if (species1 == species2) {
                        it2_start = std::next(it1);
                    }
                    for (auto it2 = it2_start; it2 != sep_list2.end(); ++it2) {
                        result.emplace_back(i, std::array{*it1, *it2});
                    }
                }
            }

            return result;
        }

        JGAP_LOG_AND_THROW("Should be unreachable");
    }

    template<size_t N, ClusterSymmetry ClusterType>
    requires (N > 1 && N <= 3)
    std::vector<SpeciesSet<N, ClusterType>> NeighbourList::getSpeciesSets() const {

        std::vector<Species> species_present;
        for (const auto &species: atoms_by_species | std::views::keys) {
            species_present.push_back(species);
        }

        std::set<SpeciesSet<N, ClusterType>> result_set;

        if constexpr (N == 2 && ClusterType == HasCentralAtom) {
            for (const auto& s0 : species_present) {
                for (const auto& s1 : species_present) {
                    result_set.emplace(s0, s1);
                }
            }
        } else if constexpr (N == 2 && ClusterType == Symmetric) {
            for (const auto& s0 : species_present) {
                for (const auto& s1 : species_present) {
                    result_set.emplace(s0, s1);
                }
            }
        } else if constexpr (N == 3 && ClusterType == HasCentralAtom) {
            for (const auto& s0 : species_present) {
                for (const auto& s1 : species_present) {
                    for (const auto& s2 : species_present) {
                        result_set.emplace(s0, s1, s2);
                    }
                }
            }
        } else if constexpr (N == 3 && ClusterType == Symmetric) {
            for (const auto& s0 : species_present) {
                for (const auto& s1 : species_present) {
                    for (const auto& s2 : species_present) {
                        result_set.emplace(s0, s1, s2);
                    }
                }
            }
        }

        return {result_set.begin(), result_set.end()};
    }

    template<size_t N>
    requires(N > 1 && N <= 3)
    std::vector<SpeciesSet<N, HasCentralAtom>> NeighbourList::getSpeciesSets(Species central_atom_species)
        const {

        if constexpr (N > 3) {
            JGAP_LOG_AND_THROW("SpeciesSet generation for N > 3 is not implemented");
        }

        std::vector<Species> species_present;
        for (const auto &species: atoms_by_species | std::views::keys) {
            species_present.push_back(species);
        }

        if (std::ranges::find(species_present, central_atom_species) == species_present.end()) {
            return {};
        }

        std::set<SpeciesSet<N, HasCentralAtom>> result_set;

        if constexpr (N == 2) {
            for (const auto& s1 : species_present) {
                result_set.emplace(central_atom_species, s1);
            }
        } else if constexpr (N == 3) {
            for (const auto& s1 : species_present) {
                for (const auto& s2 : species_present) {
                    result_set.emplace(central_atom_species, s1, s2);
                }
            }
        }

        return {result_set.begin(), result_set.end()};
    }

    template std::vector<Cluster<2, WithDerivatives>> NeighbourList
        ::findAllClusters<WithDerivatives, 2, HasCentralAtom>(const SpeciesSet<2, HasCentralAtom>& species_set) const;
    template std::vector<Cluster<2, WithDerivatives>> NeighbourList
        ::findAllClusters<WithDerivatives, 2, Symmetric>(const SpeciesSet<2, Symmetric>& species_set) const;
    template std::vector<Cluster<3, WithDerivatives>> NeighbourList
        ::findAllClusters<WithDerivatives, 3, HasCentralAtom>(const SpeciesSet<3, HasCentralAtom>& species_set) const;

    template std::vector<Cluster<2>> NeighbourList::findAllClusters<ValueOnly, 2, HasCentralAtom>(
        const SpeciesSet<2, HasCentralAtom>& species_set) const;
    template std::vector<Cluster<2>> NeighbourList::findAllClusters<ValueOnly, 2, Symmetric>(
        const SpeciesSet<2, Symmetric>& species_set) const;
    template std::vector<Cluster<3>> NeighbourList::findAllClusters<ValueOnly, 3, HasCentralAtom>(
        const SpeciesSet<3, HasCentralAtom>& species_set) const;

    template std::vector<SpeciesSet<2, HasCentralAtom>> NeighbourList::getSpeciesSets<2, HasCentralAtom>() const;
    template std::vector<SpeciesSet<2, Symmetric>> NeighbourList::getSpeciesSets<2, Symmetric>() const;
    template std::vector<SpeciesSet<3, HasCentralAtom>> NeighbourList::getSpeciesSets<3, HasCentralAtom>() const;
    template std::vector<SpeciesSet<3, Symmetric>> NeighbourList::getSpeciesSets<3, Symmetric>() const;

    template std::vector<SpeciesSet<2, HasCentralAtom>> NeighbourList::getSpeciesSets<2>(
        Species central_atom_species) const;
    template std::vector<SpeciesSet<3, HasCentralAtom>> NeighbourList::getSpeciesSets<3>(
        Species central_atom_species) const;
}