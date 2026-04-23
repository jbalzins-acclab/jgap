#include "core/fit/IsolatedAtomFit.hpp"

#include <utils/Utils.hpp>

#include "data/Vector3.hpp"
#include "core/potentials/IsolatedAtomPotential.hpp"
#include "io/log/CurrentLogger.hpp"

namespace jgap {

    std::shared_ptr<IsolatedAtomFit> IsolatedAtomFit::fromDataNode(const DataNode& params) {
        return std::make_shared<IsolatedAtomFit>(params.getOrDefault("error_on_unknown", false));
    }

    std::shared_ptr<Potential> IsolatedAtomFit::fit(const std::vector<AtomicStructure> &trainingData) {

        std::map<Species, double> isolatedEnergies = {};
        for (auto &structure: trainingData) {
            if (getOrDefault(structure.properties, "config_type", "-") == "isolated_atom") {
                if (structure.size() != 1) {
                    JGAP_LOG_AND_THROW(
                        "Structure labeled as isolated_atom does not contain exactly one atom"
                        );
                }
                if (!structure.energy.has_value()) {
                    JGAP_LOG_AND_THROW("isolated_atom with no energy");
                }

                if (const Species species = structure.species[0]; isolatedEnergies.contains(species)) {
                    if (isolatedEnergies[species] - structure.energy.value() > 1e-9) {
                        JGAP_LOG_AND_THROW(
                            "Found multiple {} isolated_atom structures with non-matching energies", species
                            );
                    }
                } else {
                    isolatedEnergies[species] = structure.energy.value();
                }
            }
        }

        return std::make_shared<IsolatedAtomPotential>(isolatedEnergies, error_on_unknown_species_);
    }
}
