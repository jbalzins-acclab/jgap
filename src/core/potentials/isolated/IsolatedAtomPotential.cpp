#include "IsolatedAtomPotential.hpp"

#include "core/atomic/Atoms.hpp"
#include "io/log/StdoutLogger.hpp"

namespace jgap {
    IsolatedAtomPotential::IsolatedAtomPotential(const std::map<Species, Real> &isolated_atom_energies)
        : isolated_energies(isolated_atom_energies) {
    }

    IsolatedAtomPotential::IsolatedAtomPotential(const std::vector<Atoms> &training_data) {

        for (const auto &atoms: training_data) {

            if (atoms.lookupConfigType().value_or("") != "isolated_atom") continue;

            if (atoms.nAtoms() != 1) {
                JGAP_LOG_AND_THROW("isolated_atom box doesn't contain one atom exactly");
            }

            if (!atoms.lookupEnergy().has_value()) {
                JGAP_LOG_AND_THROW("isolated_atom box with no energy set");
            }

            auto species = atoms.lookupSpecies()[0];

            if (isolated_energies.contains(species)) {
                JGAP_LOG_WARN("Duplicate isolated_atom box of species {}", species.symbol());
                if (abs(isolated_energies[species] - atoms.lookupEnergy().value()) > 1e-6) {
                    JGAP_LOG_AND_THROW("Energies in duplicate isolated_atom boxes of species {} vary significantly",
                                       species.symbol());
                }
            }

            isolated_energies[species] = atoms.lookupEnergy().value();
        }
    }

    AtomicQuantity IsolatedAtomPotential::calculateEnergy(const Atoms &atoms) const {

        AtomicQuantity result(atoms.nAtoms());

        for (const auto &species: atoms.lookupSpecies()) {
            if (isolated_energies.contains(species)) {
                result.value += isolated_energies.at(species);
            }
        }

        return result;
    }
}
