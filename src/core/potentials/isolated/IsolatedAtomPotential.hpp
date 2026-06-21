#ifndef JGAP_ISOLATEDATOMPOTENTIAL_HPP
#define JGAP_ISOLATEDATOMPOTENTIAL_HPP

#include <map>
#include "../Potential.hpp"
#include "core/atomic/species/Species.hpp"
#include "../../../serialization/SerializationRegistry.hpp"

namespace jgap {
    class IsolatedAtomPotential : public Potential {
    public:
        IsolatedAtomPotential(const std::map<Species, Real>& isolated_atom_energies);
        IsolatedAtomPotential(const std::vector<Atoms>& training_data);

        AtomicQuantity calculateEnergy(const Atoms &atoms) const override;

        Cutoffs getCutoffs() const override { return {}; }

        void fillTables(TabulationData &table) const override {
            for (auto& [species, energy]: isolated_energies) {
                table.isolated_energies[species] += energy;
            }
        }

        IsolatedAtomPotential* clone() const override {
            return new IsolatedAtomPotential(*this);
        }

        const std::map<Species, Real>& getIsolatedEnergies() const {
            return isolated_energies;
        }

    private:
        std::map<Species, Real> isolated_energies;
    };
}

#endif