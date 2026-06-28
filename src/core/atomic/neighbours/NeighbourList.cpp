#include "NeighbourList.hpp"
#include <tbb/parallel_for.h>
#include <iostream>
#include <algorithm>
#include <numeric>
#include <ranges>
#include <set>
#include "io/log/CurrentLogger.hpp"

namespace jgap {

    std::array<int, 3> NeighbourList::findMaxRep(const Atoms& structure, const Real cutoff) {
        auto pbc = structure.lookupPbc();
        auto lattice = structure.lookupLattice();

        if (!lattice.has_value()) {
            if (pbc[0] || pbc[1] || pbc[2]) {
                JGAP_LOG_AND_THROW("PBC is true but no lattice is provided.");
            }
            return {0, 0, 0};
        }

        std::array maxRep = {
            pbc[0] ? static_cast<int>(cutoff / lattice->a.norm() + 1) : 0,
            pbc[1] ? static_cast<int>(cutoff / lattice->b.norm() + 1) : 0,
            pbc[2] ? static_cast<int>(cutoff / lattice->c.norm() + 1) : 0
        };

        // triclinic
        if (abs(lattice->a.dot(lattice->b)) > 1e-6
            || abs(lattice->a.dot(lattice->c)) > 1e-6
            || abs(lattice->b.dot(lattice->c)) > 1e-6)
        {
            if (pbc[0]) maxRep[0] = static_cast<int>(cutoff / lattice->a.aproject(lattice->b, lattice->c)) + 1;
            if (pbc[1]) maxRep[1] = static_cast<int>(cutoff / lattice->b.aproject(lattice->a, lattice->c)) + 1;
            if (pbc[2]) maxRep[2] = static_cast<int>(cutoff / lattice->c.aproject(lattice->b, lattice->a)) + 1;
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

        Vector3 a{}, b{}, c{};
        if (lattice_opt.has_value()) {
            a = lattice_opt->a;
            b = lattice_opt->b;
            c = lattice_opt->c;
        }

        for (size_t i = 0; i < box.nAtoms(); i++) {
            auto& neighbours_i = neighbours_per_atom[i];
            auto pos_i = positions[i];

            for (size_t j = 0; j < box.nAtoms(); j++) {
                auto pos_j = positions[j];

                for (int rep_a = -max_rep[0]; rep_a <= max_rep[0]; rep_a++) {
                    for (int rep_b = -max_rep[1]; rep_b <= max_rep[1]; rep_b++) {
                        for (int rep_c = -max_rep[2]; rep_c <= max_rep[2]; rep_c++) {
                            if (i == j && rep_a == 0 && rep_b == 0 && rep_c == 0) continue;

                            Vector3 offset = a * rep_a + b * rep_b + c * rep_c;

                            auto separation = Separation(pos_i, pos_j + offset);

                            if (separation.magnitude <= cutoff) {
                                neighbours_i[species[j]].emplace_back(j, separation);
                            }
                        }
                    }
                }
            }
        }
    }

    std::vector<NeighbourList> NeighbourList::form(const std::vector<Atoms>& boxes, Real cutoff) {
        std::vector<NeighbourList> result;
        result.reserve(boxes.size());

        tbb::parallel_for(static_cast<size_t>(0), boxes.size(), [&](size_t i) {
            result[i] = NeighbourList(boxes[i], cutoff);
        });

        return result;
    }

    // CalcType first since template deduction works well for others
    /*template<CalculationType CalcType, size_t N, ClusterSymmetry ClusterType>
    requires(N > 1 && N <= 3)
    std::vector<Cluster<N, CalcType>> NeighbourList::findAllClusters(
        const SpeciesSet<N, ClusterType> &species_set) const {

    }*/


    template<size_t N, ClusterSymmetry ClusterType>
    requires (N > 1 && N <= 3)
    std::set<SpeciesSet<N, ClusterType>> NeighbourList::getSpeciesSets() const {

        std::vector<Species> species_present;
        for (const auto &species: atoms_by_species | std::views::keys) {
            species_present.push_back(species);
        }

        std::set<SpeciesSet<N, ClusterType>> result_set;

        if constexpr (N == 2 && ClusterType == NodeSymmetric) {
            for (const auto& s0: species_present) {
                for (const auto& s1: species_present) {
                    result_set.emplace(s0, s1);
                }
            }
        } else if constexpr (N == 2 && ClusterType == FullSymmetry) {
            for (const auto& s0: species_present) {
                for (const auto& s1: species_present) {
                    result_set.emplace(s0, s1);
                }
            }
        } else if constexpr (N == 3 && ClusterType == NodeSymmetric) {
            for (const auto& s0: species_present) {
                for (const auto& s1: species_present) {
                    for (const auto& s2: species_present) {
                        result_set.emplace(s0, s1, s2);
                    }
                }
            }
        } else if constexpr (N == 3 && ClusterType == FullSymmetry) {
            for (const auto& s0: species_present) {
                for (const auto& s1: species_present) {
                    for (const auto& s2: species_present) {
                        result_set.emplace(s0, s1, s2);
                    }
                }
            }
        }

        return result_set;
    }

    template<size_t N>
    requires(N > 1 && N <= 3)
    std::set<SpeciesSet<N, NodeSymmetric>> NeighbourList::getSpeciesSets(Species central_atom_species) const {

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

        std::set<SpeciesSet<N, NodeSymmetric>> result_set;

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

        return result_set;
    }

    template std::vector<Cluster<2, WithGradients>> NeighbourList
        ::findAllClusters<WithGradients, 2, NodeSymmetric>(const SpeciesSet<2, NodeSymmetric>& species_set) const;
    template std::vector<Cluster<2, WithGradients>> NeighbourList
        ::findAllClusters<WithGradients, 2, FullSymmetry>(const SpeciesSet<2, FullSymmetry>& species_set) const;
    template std::vector<Cluster<3, WithGradients>> NeighbourList
        ::findAllClusters<WithGradients, 3, NodeSymmetric>(const SpeciesSet<3, NodeSymmetric>& species_set) const;

    template std::vector<Cluster<2>> NeighbourList::findAllClusters<ValueOnly, 2, NodeSymmetric>(
        const SpeciesSet<2, NodeSymmetric>& species_set) const;
    template std::vector<Cluster<2>> NeighbourList::findAllClusters<ValueOnly, 2, FullSymmetry>(
        const SpeciesSet<2, FullSymmetry>& species_set) const;
    template std::vector<Cluster<3>> NeighbourList::findAllClusters<ValueOnly, 3, NodeSymmetric>(
        const SpeciesSet<3, NodeSymmetric>& species_set) const;

    template std::set<SpeciesSet<2, NodeSymmetric>> NeighbourList::getSpeciesSets<2, NodeSymmetric>() const;
    template std::set<SpeciesSet<2, FullSymmetry>> NeighbourList::getSpeciesSets<2, FullSymmetry>() const;
    template std::set<SpeciesSet<3, NodeSymmetric>> NeighbourList::getSpeciesSets<3, NodeSymmetric>() const;
    template std::set<SpeciesSet<3, FullSymmetry>> NeighbourList::getSpeciesSets<3, FullSymmetry>() const;

    template std::set<SpeciesSet<2, NodeSymmetric>> NeighbourList::getSpeciesSets<2>(Species central_atom_species)
        const;
    template std::set<SpeciesSet<3, NodeSymmetric>> NeighbourList::getSpeciesSets<3>(Species central_atom_species)
        const;
}