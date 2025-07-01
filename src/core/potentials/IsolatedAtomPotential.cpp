//
// Created by Jegors Balzins on 22.6.2025.
//

#include "core/potentials/IsolatedAtomPotential.hpp"

#include "io/log/StdoutLogger.hpp"

namespace jgap {
    IsolatedAtomPotential::IsolatedAtomPotential(const vector<AtomicStructure> &containsIsolatedAtoms,
                                                 const bool errorOnUnknownSpecies)
        : _errorOnUnknownSpecies(errorOnUnknownSpecies) {

        StdoutLogger::initIfNotInitialized();

        _isolatedEnergies = {};
        for (auto &structure: containsIsolatedAtoms) {
            if (structure.configType.value_or("-") == "isolated_atom") {
                if (structure.atoms.size() != 1) {
                    Logger::logger -> error(
                        "Structure labeled as isolated_atom does not contain exactly one atom",
                        true
                        );
                }
                if (!structure.energy.has_value()) {
                    Logger::logger -> error("isolated_atom with no energy",true);
                }

                const Species species = structure.atoms[0].species;

                if (_isolatedEnergies.contains(species)) {
                    if (_isolatedEnergies[species] - structure.energy.value() > 1e-9) {
                        Logger::logger -> error(
                            format("Found multiple {} isolated_atom structures with non-matching energies", species),
                            true
                            );
                    }
                } else {
                    _isolatedEnergies[species] = structure.energy.value();
                }
            }
        }
    }

    PotentialPrediction IsolatedAtomPotential::predict(const AtomicStructure& structure) {

        double result = 0;

        for (const auto &atom: structure.atoms) {
            const auto species = atom.species;

            if (_isolatedEnergies.contains(species)) {
                result += _isolatedEnergies[species];
            } else {
                if (_errorOnUnknownSpecies) {
                    Logger::logger -> error("Unknown isolated_atom energy for " + species,true);
                }
            }
        }

        return PotentialPrediction{.energy = result};
    }

    nlohmann::json IsolatedAtomPotential::serialize() {
        return _isolatedEnergies;
    }
}
