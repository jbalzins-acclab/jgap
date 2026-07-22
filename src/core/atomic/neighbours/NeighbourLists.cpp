#include "NeighbourLists.hpp"
#include <algorithm>
#include <iostream>
#include <numeric>
#include <ranges>
#include <set>

#include "core/UnseqFor.hpp"
#include "io/log/CurrentLogger.hpp"

namespace jgap {

    std::array<int, 3> NeighbourLists::findMaxRep(const Atoms& structure, const Real cutoff) {
        auto pbc = structure.getPbc();
        auto lattice = structure.getLattice();

        if (!lattice.has_value()) {
            if (pbc[0] || pbc[1] || pbc[2]) {
                JGAP_LOG_AND_THROW("PBC is true but no lattice is provided.");
            }
            return {0, 0, 0};
        }

        std::array maxRep = {pbc[0] ? static_cast<int>(cutoff / lattice->a.norm() + 1) : 0,
                             pbc[1] ? static_cast<int>(cutoff / lattice->b.norm() + 1) : 0,
                             pbc[2] ? static_cast<int>(cutoff / lattice->c.norm() + 1) : 0};

        // triclinic
        if (abs(lattice->a.dot(lattice->b)) > 1e-6 || abs(lattice->a.dot(lattice->c)) > 1e-6 ||
            abs(lattice->b.dot(lattice->c)) > 1e-6) {
            if (pbc[0]) maxRep[0] = static_cast<int>(cutoff / lattice->a.aproject(lattice->b, lattice->c)) + 1;
            if (pbc[1]) maxRep[1] = static_cast<int>(cutoff / lattice->b.aproject(lattice->a, lattice->c)) + 1;
            if (pbc[2]) maxRep[2] = static_cast<int>(cutoff / lattice->c.aproject(lattice->b, lattice->a)) + 1;
        }

        return maxRep;
    }

    NeighbourLists::NeighbourLists(const Atoms& box, Real cutoff) : cutoff(cutoff) {
        const auto max_rep = findMaxRep(box, cutoff);
        auto lattice_opt = box.getLattice();
        auto& species = box.getSpecies();
        auto& positions = box.getPositions();

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
                            if (i == j && rep_a == 0 && rep_b == 0 && rep_c == 0) [[unlikely]]
                                continue;

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

    std::vector<NeighbourLists> NeighbourLists::form(const std::vector<Atoms>& boxes, Real cutoff) {
        std::vector<std::optional<NeighbourLists>> result_opts(boxes.size());

        unseqForIndex(0, boxes.size(),
                      [&boxes, &result_opts, cutoff](size_t i) { result_opts[i] = NeighbourLists(boxes[i], cutoff); });

        return std::views::transform(result_opts,
                                     [](std::optional<NeighbourLists>& opt) { return std::move(opt.value()); }) |
               std::ranges::to<std::vector>();
    }
}
