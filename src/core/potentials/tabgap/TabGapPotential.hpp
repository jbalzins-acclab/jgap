#ifndef JGAP_TABGAPPOTENTIAL_HPP
#define JGAP_TABGAPPOTENTIAL_HPP

#include "components/EamTGComponent.hpp"
#include "components/TabGapComponent.hpp"
#include "core/potentials/Potential.hpp"

namespace jgap {

    class TabGapPotential : public Potential {
    public:
        friend class TabGapIO;

        TabGapPotential(TabulationData energy_tables);

        AtomicQuantity calculateEnergy(const Atoms& atoms) const override;

        Cutoffs getCutoffs() const override;

        void fillTables(TabulationData& table) const override;

        std::unique_ptr<Potential> clone() const override {
            return std::make_unique<TabGapPotential>(*this);
        }

    private:
        size_t n_2b_components{};
        size_t n_3b_components{};
        size_t n_eam_components{};

        std::map<Species, Real> isolated_atom_energies;
        std::vector<ValuePtr<TabGapComponent>> components;

        TabGapPotential(const std::map<Species, Real>& isolated_atom_energies = {},
                        const std::vector<ValuePtr<TabGapComponent>>& components = {});

        // Recomputes n_2b/n_3b/n_eam from the current components (the reading paths build `components`
        // directly, so the counts that the writer relies on must be derived afterwards).
        void recomputeComponentCounts();
    };
}

#endif
