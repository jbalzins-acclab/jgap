#include "GapPotential.hpp"

#include <numeric>

namespace jgap {

    void GapPotential::setCoefficients(const std::vector<Real> &new_coefficients) {
        if (new_coefficients.empty()) {
            JGAP_LOG_AND_THROW("Cannot set empty coefficients");
        }

        size_t n = std::transform_reduce(
            components.begin(),
            components.end(),
            size_t{0},
            std::plus(),
            [](const auto& ptr) { return ptr->nSparsePoints(); }
        );

        if (n != new_coefficients.size()) {
            JGAP_LOG_AND_THROW("{} sparse points, but {} coefficients", n, new_coefficients.size());
        }

        auto iter = new_coefficients.begin();
        for (auto& component: components) {
            component->setCoefficients(iter);
        }
    }

    AtomicQuantity GapPotential::calculateEnergy(const Atoms &atoms) const {
        AtomicQuantity result(atoms.nAtoms());

        const NeighbourList neighbour_list(atoms, getCutoffs().maxOverall());

        if (optional_external_potential.get() != nullptr) {
            result += optional_external_potential->calculateEnergy(atoms);
        }

        for (const auto& component: components) {
            result += component->energy(neighbour_list);
        }

        return result;
    }

    Cutoffs GapPotential::getCutoffs() const {
        Cutoffs result{};

        if (optional_external_potential.get() != nullptr) {
            result = optional_external_potential->getCutoffs();
        }

        for (const auto& component: components) {
            result += component->getCutoffs();
        }

        return result;
    }

    const std::vector<ValuePtr<GapComponent>>& GapPotential::getComponents() const {
        return components;
    }

    void GapPotential::fillTables(TabulationData &table) const {
        JGAP_LOG_INFO("Tabulating GAP");

        JGAP_LOG_DEBUG("Tabulating external potential");
        if (optional_external_potential.get() != nullptr) {
            optional_external_potential->fillTables(table);
        }

        JGAP_LOG_DEBUG("Tabulating components");
        for (const auto& component: components) {
            component->tabulate(table);
        }

        JGAP_LOG_INFO("Finished tabulating GAP");
    }
}
