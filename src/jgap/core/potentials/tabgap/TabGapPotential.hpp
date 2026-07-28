#ifndef JGAP_TABGAPPOTENTIAL_HPP
#define JGAP_TABGAPPOTENTIAL_HPP

#include <map>
#include "components/EamTGComponent.hpp"
#include "components/ThreeBodyTGComponent.hpp"
#include "components/TwoBodyTGComponent.hpp"
#include "jgap/core/potentials/Potential.hpp"

namespace jgap {

    class TabGapPotential : public Potential {
    public:
        TabGapPotential() = default;
        explicit TabGapPotential(TabulationData energy_tables);

        TabGapPotential(
            std::map<Species, Real> isolated_atom_energies,
            std::map<Species2Sorted, TwoBodyTGComponent> two_body_components = {},
            std::map<Species3AtomicSorted, ThreeBodyTGComponent> three_body_components = {},
            std::multimap<Species, EamTGComponent> eam_components = {}
        );

        AtomicQuantity calculateEnergy(const Atoms& atoms) const override;

        Cutoffs getCutoffs() const override;

        void fillTables(TabulationData& table) const override;

        TabGapPotential* clone() const override { return new TabGapPotential(*this); }

        const auto& getIsolatedAtomEnergies() const { return isolated_atom_energies; }
        const auto& getTwoBodyComponents() const { return two_body_components; }
        const auto& getThreeBodyComponents() const { return three_body_components; }
        const auto& getEamComponents() const { return eam_components; }

    private:
        std::map<Species, Real> isolated_atom_energies;
        std::map<Species2Sorted, TwoBodyTGComponent> two_body_components;
        std::map<Species3AtomicSorted, ThreeBodyTGComponent> three_body_components;
        std::multimap<Species, EamTGComponent> eam_components;
    };
}

#endif
