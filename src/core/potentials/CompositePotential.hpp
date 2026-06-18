#ifndef JGAP_COMPOSITEPOTENTIAL_HPP
#define JGAP_COMPOSITEPOTENTIAL_HPP

#include "Potential.hpp"
#include <map>
#include <vector>
#include <memory>
#include <string>
#include <utility>

#include "../ValuePtr.hpp"

namespace jgap {

    class CompositePotential final : public Potential {
    public:
        CompositePotential(std::map<std::string, ValuePtr<Potential>> potentials_map);
        CompositePotential(std::vector<ValuePtr<Potential>> potentials_list);

        template<typename... Potentials>
        CompositePotential(Potentials... pots) {
            int i = 0;
            ( (potentials[std::to_string(i++)] = std::move(pots)), ... );
        }

        Cutoffs getCutoffs() const override;

        AtomicQuantity calculateEnergy(const Atoms &atoms) const override;

        const std::map<std::string, ValuePtr<Potential>>& getPotentials() const { return potentials; }
        std::map<std::string, ValuePtr<Potential>>& getPotentials() { return potentials; }

        void fillTables(TabulationData &table) const override;

        std::unique_ptr<Potential> clone() const override {
            return std::make_unique<CompositePotential>(*this);
        }

    private:
        std::map<std::string, ValuePtr<Potential>> potentials;
    };
}

#endif