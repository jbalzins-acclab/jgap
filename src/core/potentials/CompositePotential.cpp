#include "core/potentials/CompositePotential.hpp"
#include <ranges>

#include "core/atomic/Atoms.hpp"

namespace jgap {

    CompositePotential::CompositePotential(std::map<std::string, std::unique_ptr<Potential>> potentials_map)
        : potentials(std::move(potentials_map)) {
    }

    CompositePotential::CompositePotential(std::vector<std::unique_ptr<Potential>> potentials_list) {
        for (size_t i = 0; i < potentials_list.size(); ++i) {
            potentials[std::to_string(i + 1)] = std::move(potentials_list[i]);
        }
    }

    Cutoffs CompositePotential::getCutoffs() {
        Cutoffs res{};
        for (const auto &potential: potentials | std::views::values) {
            res += potential->getCutoffs();
        }
        return res;
    }

    AtomicQuantity CompositePotential::calculateEnergy(const Atoms &atoms) {
        AtomicQuantity result(atoms.nAtoms());
        for (const auto &potential : potentials | std::views::values) {
            result += potential->calculateEnergy(atoms);
        }
        return result;
    }

}
