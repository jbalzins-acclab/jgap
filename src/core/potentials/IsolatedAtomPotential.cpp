//
// Created by Jegors Balzins on 22.6.2025.
//

#include "core/potentials/IsolatedAtomPotential.hpp"

#include "io/log/StdoutLogger.hpp"

namespace jgap {
    IsolatedAtomPotential::IsolatedAtomPotential(const nlohmann::json& params) {
        _errorOnUnknownSpecies = params.value("error_on_unknown", true);
        _isolatedEnergies = {};
        for (const auto& [element, energy]: params["energies"].items()) {
            _isolatedEnergies[element] = energy.get<double>();
        }
    }

    nlohmann::json IsolatedAtomPotential::serialize() {
        return {
            {"error_on_unknown", _errorOnUnknownSpecies},
            {"energies", _isolatedEnergies}
        };
    }

    PotentialPrediction IsolatedAtomPotential::predict(const AtomicStructure& structure) {

        double result = 0;

        for (const auto &atom: structure.atoms) {
            const auto species = atom.species;

            if (_isolatedEnergies.contains(species)) {
                result += _isolatedEnergies[species];
            } else {
                if (_errorOnUnknownSpecies) {
                    CurrentLogger::get() -> error("Unknown isolated_atom energy for " + species,true);
                }
            }
        }

        return PotentialPrediction{.energy = result};
    }
}
