#ifndef ISOLATEDATOMPOTENTIAL_HPP
#define ISOLATEDATOMPOTENTIAL_HPP

#include <map>
#include "../Potential.hpp"
#include "core/atomic/species/Species.hpp"
#include "io/Serializable.hpp"
#include "io/parse/ParserRegistry.hpp"

namespace jgap {
    class IsolatedAtomPotential : public Potential {
    public:
        IsolatedAtomPotential(const std::map<Species, double>& isolated_atom_energies, bool error_on_unknown);
        IsolatedAtomPotential(const std::vector<Atoms>& training_data, bool error_on_unknown);

        AtomicQuantity calculateEnergy(const Atoms &atoms) override;

        Cutoffs getCutoffs() override { return {}; }

    private:
        bool error_on_unknown_species;
        std::map<Species, double> isolated_energies;
    };
}

#endif
