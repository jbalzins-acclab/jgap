#ifndef ISOLATEDATOMPOTENTIAL_HPP
#define ISOLATEDATOMPOTENTIAL_HPP

#include <map>
#include "../Potential.hpp"
#include "core/atomic/species/Species.hpp"
#include "io/parse/ParserRegistry.hpp"

namespace jgap {
    class IsolatedAtomPotential : public Potential {
    public:
        IsolatedAtomPotential(const std::map<Species, Real>& isolated_atom_energies, bool error_on_unknown);
        IsolatedAtomPotential(const std::vector<Atoms>& training_data, bool error_on_unknown);

        AtomicQuantity calculateEnergy(const Atoms &atoms) const override;

        Cutoffs getCutoffs() const override { return {}; }

        void tabulate(TabulationData &table) const override {
            for (auto& [species, energy]: isolated_energies) {
                table.isolated_energies[species] += energy;
            }
        }

        std::unique_ptr<Potential> clone() const override {
            return std::make_unique<IsolatedAtomPotential>(*this);
        }

    private:
        bool error_on_unknown_species;
        std::map<Species, Real> isolated_energies;
    };
}

#endif
