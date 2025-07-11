//
// Created by Jegors Balzins on 9.7.2025.
//

#include "core/fit/IsolatedAtomFit.hpp"

#include "core/potentials/IsolatedAtomPotential.hpp"

namespace jgap {

    IsolatedAtomFit::IsolatedAtomFit(const nlohmann::json& params) {
        _errorOnUnknownSpecies = params.value("error_on_unknown", true);
    }

    shared_ptr<Potential> IsolatedAtomFit::fit(const vector<AtomicStructure> &trainingData) {

        map<Species, double> isolatedEnergies = {};
        for (auto &structure: trainingData) {
            if (structure.configType.value_or("-") == "isolated_atom") {
                if (structure.atoms.size() != 1) {
                    CurrentLogger::get() -> error(
                        "Structure labeled as isolated_atom does not contain exactly one atom",
                        true
                        );
                }
                if (!structure.energy.has_value()) {
                    CurrentLogger::get() -> error("isolated_atom with no energy",true);
                }

                const Species species = structure.atoms[0].species;

                if (isolatedEnergies.contains(species)) {
                    if (isolatedEnergies[species] - structure.energy.value() > 1e-9) {
                        CurrentLogger::get() -> error(
                            format("Found multiple {} isolated_atom structures with non-matching energies", species),
                            true
                            );
                    }
                } else {
                    isolatedEnergies[species] = structure.energy.value();
                }
            }
        }

        return make_shared<IsolatedAtomPotential>(nlohmann::json{
            {"error_on_unknown", _errorOnUnknownSpecies},
            {"energies", isolatedEnergies}
        });
    }
}
