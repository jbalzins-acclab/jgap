#include "core/potentials/CompositePotential.hpp"
#include <ranges>

#include "core/atomic/Atoms.hpp"

namespace jgap {

    CompositePotential::CompositePotential(std::map<std::string, ValuePtr<Potential>> potentials_map)
        : potentials(std::move(potentials_map)) {
    }

    CompositePotential::CompositePotential(std::vector<ValuePtr<Potential>> potentials_list) {
        for (size_t i = 0; i < potentials_list.size(); ++i) {
            potentials[std::to_string(i + 1)] = std::move(potentials_list[i]);
        }
    }

    Cutoffs CompositePotential::getCutoffs() const {
        Cutoffs res{};
        for (const auto &potential: potentials | std::views::values) {
            res += potential->getCutoffs();
        }
        return res;
    }

    AtomicQuantity CompositePotential::calculateEnergy(const Atoms &atoms) const {
        AtomicQuantity result(atoms.nAtoms());
        for (const auto &potential : potentials | std::views::values) {
            result += potential->calculateEnergy(atoms);
        }
        return result;
    }

    void CompositePotential::fillTables(TabulationData &table) const {
        for (const auto &potential : potentials | std::views::values) {
            potential->fillTables(table);
        }
    }
}
