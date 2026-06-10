#ifndef JGAP_COMPOSITEPOTENTIAL_HPP
#define JGAP_COMPOSITEPOTENTIAL_HPP

#include "Potential.hpp"
#include <map>
#include <vector>
#include <memory>
#include <string>
#include <utility>

namespace jgap {

    class CompositePotential final : public Potential {
    public:
        CompositePotential(std::map<std::string, std::unique_ptr<Potential>> potentials_map);
        CompositePotential(std::vector<std::unique_ptr<Potential>> potentials_list);

        template<typename... Potentials>
        CompositePotential(Potentials... pots) {
            int i = 0;
            ( (potentials[std::to_string(i++)] = std::move(pots)), ... );
        }

        Cutoffs getCutoffs() override;
        AtomicQuantity calculateEnergy(const Atoms &atoms) override;

        const std::map<std::string, std::unique_ptr<Potential>>& getPotentials() const { return potentials; }
        std::map<std::string, std::unique_ptr<Potential>>& getPotentials() { return potentials; }

    private:
        std::map<std::string, std::unique_ptr<Potential>> potentials;
    };
}

#endif