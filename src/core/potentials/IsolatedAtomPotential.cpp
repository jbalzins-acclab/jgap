//
// Created by Jegors Balzins on 22.6.2025.
//

#include "core/potentials/IsolatedAtomPotential.hpp"

#include "io/log/StdoutLogger.hpp"

namespace jgap {
    IsolatedAtomPotential::IsolatedAtomPotential(const nlohmann::json& params) {
        CurrentLogger::get()->debug("Parsing isolated atom potentials params");
        _errorOnUnknownSpecies = params.value("error_on_unknown", true);
        _isolatedEnergies = {};
        for (const auto& [element, energy]: params["energies"].items()) {
            _isolatedEnergies[element] = energy.get<double>();
        }
    }

    IsolatedAtomPotential::IsolatedAtomPotential(const map<Species, double> &isolatedAtomEnergies, bool errorOnUnknown) {
        _isolatedEnergies = isolatedAtomEnergies;
        _errorOnUnknownSpecies = errorOnUnknown;
    }

    nlohmann::json IsolatedAtomPotential::serialize() {
        return {
            {"error_on_unknown", _errorOnUnknownSpecies},
            {"energies", _isolatedEnergies}
        };
    }

    PotentialPrediction IsolatedAtomPotential::predict(const AtomicStructure& structure) {

        double result = 0;

        for (const auto &atom: structure) {
            if (_isolatedEnergies.contains(atom.species())) {
                result += _isolatedEnergies[atom.species()];
            } else {
                if (_errorOnUnknownSpecies) {
                    CurrentLogger::get() -> error("Unknown isolated_atom energy for " + atom.species(),true);
                }
            }
        }

        return PotentialPrediction{.energy = result};
    }

    TabulationData IsolatedAtomPotential::tabulate(const TabulationParams &params) {
        TabulationData result;
        result.isolatedEnergies = _isolatedEnergies;
        return result;
    }
}
