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
        for (const auto& [species, atom_indices] : atoms_by_species) {
            for (size_t index : atom_indices) {
                if (index == atom_index) {
                    return species;
                }
            }
        }
        throw std::runtime_error("Atom index not found in any species.");
    }

    std::array<int, 3> NeighbourList::findMaxRep(const Atoms& structure, const Real cutoff) {
        auto pbc = structure.getPbc();
        auto lattice_opt = structure.getLattice();

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
            pbc[0] ? static_cast<int>(cutoff / side1.norm() + 2) : 0,
            pbc[1] ? static_cast<int>(cutoff / side2.norm() + 2) : 0,
            pbc[2] ? static_cast<int>(cutoff / side3.norm() + 2) : 0
        };

        // triclinic
        if (abs(side1.dot(side2)) > 1e-6 || abs(side1.dot(side3)) > 1e-6 || abs(side2.dot(side3)) > 1e-6) {
            if (pbc[0]) maxRep[0] = static_cast<int>(cutoff / side1.aproject(side2, side3)) + 2;
            if (pbc[1]) maxRep[1] = static_cast<int>(cutoff / side2.aproject(side1, side3)) + 2;
            if (pbc[2]) maxRep[2] = static_cast<int>(cutoff / side3.aproject(side2, side1)) + 2;
        }

        return maxRep;
    }

    NeighbourList::NeighbourList(const Atoms& box, Real cutoff) : cutoff(cutoff) {
        const auto max_rep = findMaxRep(box, cutoff);
        auto lattice_opt = box.getLattice();

        neighbours_per_atom.resize(box.nAtoms());

        for (size_t i = 0; i < box.nAtoms(); i++) {
            atoms_by_species[box.getSpecies()[i]].push_back(i);
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

                            auto pos_j_plus_offset = box.getPositions()[j] + offset;
                            auto separation_ij = Separation(box.getPositions()[i], pos_j_plus_offset);

                            if (separation_ij.magnitude > 0 && separation_ij.magnitude <= cutoff) {
                                neighbours_i[box.getSpecies()[j]].emplace_back(j, separation_ij);
                            }
                        }
                    }
                }
            }
            // Sort neighbours by distance
            for (auto &val: neighbours_i | std::views::values) {
                std::sort(val.begin(), val.end(),
                    [](const NeighbourData& a, const NeighbourData& b) {
                        return a.separation.magnitude < b.separation.magnitude;
                });
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

    template<size_t N>
    requires(N > 1)
    std::vector<Separations<N>> NeighbourList::findAllSeparations(const SpeciesSet &species_set,
                                                                  std::optional<Real> max_distance) const {

        if constexpr (N > 3) {
            JGAP_LOG_AND_THROW("Inter-separations finding for N > 3 is not implemented");
        }

        std::array<std::optional<Species>, N> required_species;
        for (size_t i = 0; i < N; ++i) {
            int16_t target_id = species_set[i];
            if (target_id == -1) continue;
            for (const auto& pair : atoms_by_species) {
                if (pair.first.id() == target_id) {
                    required_species[i] = pair.first;
                    break;
                }
            }
        }

        // For N=2 and N=3, check if we found the required species
        for (size_t i = 0; i < N; ++i) {
            if (!required_species[i].has_value()) return {};
        }

        const Real effective_cutoff = max_distance.value_or(cutoff);

        if constexpr (N == 2) {

            std::vector<Separations<2>> result;

            const auto& species0 = required_species[0].value();
            const auto& species1 = required_species[1].value();
            if (atoms_by_species.count(species0)) {
                for (size_t i : atoms_by_species.at(species0)) {
                    if (neighbours_per_atom[i].count(species1)) {
                        const auto& list = neighbours_per_atom[i].at(species1);
                        auto it_end = list.end();
                        if (max_distance.has_value()) {
                            it_end = std::partition_point(list.begin(), list.end(),
                                [&](const NeighbourData& nd) {
                                    return nd.separation.magnitude <= effective_cutoff;
                            });
                        }
                        for (auto it = list.begin(); it != it_end; ++it) {
                            Separations<2> seps;
                            seps.atom_indexes[0] = i;
                            seps.atom_indexes[1] = it->atom_index;
                            seps.separations[0] = it->separation;
                            result.push_back(seps);
                        }
                    }
                }
            }

            return result;

        } else if constexpr (N == 3) {

            std::vector<Separations<3>> result;

            const auto& species0 = required_species[0].value();
            const auto& species1 = required_species[1].value();
            const auto& species2 = required_species[2].value();

            if (atoms_by_species.count(species0)) {
                for (size_t i : atoms_by_species.at(species0)) {
                    if (neighbours_per_atom[i].count(species1) && neighbours_per_atom[i].count(species2)) {
                        const auto& list1_full = neighbours_per_atom[i].at(species1);
                        const auto& list2_full = neighbours_per_atom[i].at(species2);

                        auto it1_end = list1_full.end();
                        auto it2_end = list2_full.end();

                        if (max_distance.has_value()) {
                            it1_end = std::partition_point(list1_full.begin(), list1_full.end(),
                                [&](const NeighbourData& nd) {
                                    return nd.separation.magnitude <= effective_cutoff;
                            });
                            it2_end = std::partition_point(list2_full.begin(), list2_full.end(),
                                [&](const NeighbourData& nd) {
                                    return nd.separation.magnitude <= effective_cutoff;
                            });
                        }

                        for (auto it1 = list1_full.begin(); it1 != it1_end; ++it1) {
                            size_t idx1 = std::distance(list1_full.begin(), it1);
                            auto it2_start = list2_full.begin();
                            if (species1 == species2) {
                                it2_start = it1 + 1;
                            }
                            for (auto it2 = it2_start; it2 != it2_end; ++it2) {
                                Separations<3> seps;
                                seps.atom_indexes[0] = i;
                                seps.atom_indexes[1] = it1->atom_index;
                                seps.atom_indexes[2] = it2->atom_index;
                                seps.separations[0] = it1->separation; // 0-1
                                seps.separations[1] = it2->separation; // 0-2
                                seps.separations[2] = Separation(
                                    it1->separation.vec(),
                                    it2->separation.vec()
                                    ); // 1-2
                                result.push_back(seps);
                            }
                        }
                    }
                }
            }

            return result;
        }
        return {};
    }

    template<size_t N>
    requires(N > 0)
    std::vector<SpeciesSet> NeighbourList::getSpeciesSets() const {
        if constexpr (N > 3) {
            JGAP_LOG_AND_THROW("SpeciesSet generation for N > 3 is not implemented");
        }

        std::set<Species> species_present;
        for (const auto &species: atoms_by_species | std::views::keys) {
            species_present.insert(species);
        }

        std::set<SpeciesSet> result_set;

        if constexpr (N == 1) {
            for (const auto& s : species_present) {
                result_set.insert(SpeciesSet(s));
            }
        } else if constexpr (N == 2) {
            for (const auto& s0 : species_present) {
                for (const auto& s1 : species_present) {
                    SpeciesSet ss(s0);
                    ss.node(s1);
                    result_set.insert(ss);
                }
            }
        } else if constexpr (N == 3) {
            for (const auto& s0 : species_present) {
                for (const auto& s1 : species_present) {
                    for (const auto& s2 : species_present) {
                        SpeciesSet ss(s0);
                        ss.node(s1).node(s2);
                        result_set.insert(ss);
                    }
                }
            }
        }

        return {result_set.begin(), result_set.end()};
    }

    template std::vector<Separations<2>> NeighbourList::findAllSeparations(const SpeciesSet& species_set,
                                                                           std::optional<Real> max_distance) const;
    template std::vector<Separations<3>> NeighbourList::findAllSeparations(const SpeciesSet& species_set,
                                                                           std::optional<Real> max_distance) const;

    template std::vector<SpeciesSet> NeighbourList::getSpeciesSets<1>() const;
    template std::vector<SpeciesSet> NeighbourList::getSpeciesSets<2>() const;
    template std::vector<SpeciesSet> NeighbourList::getSpeciesSets<3>() const;

}
