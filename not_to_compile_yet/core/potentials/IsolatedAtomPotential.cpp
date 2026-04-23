#include "core/potentials/IsolatedAtomPotential.hpp"

#include "io/log/StdoutLogger.hpp"

namespace jgap {
    IsolatedAtomPotential::IsolatedAtomPotential(const DataNode& params) {
        JGAP_LOG_DEBUG("Parsing isolated atom potentials params");
        error_on_unknown_species_ = params.value("error_on_unknown", true);
        isolated_energies_ = {};
        for (const auto& [element, energy]: params["energies"].items()) {
            isolated_energies_[element] = energy.get<double>();
        }
    }

    IsolatedAtomPotential::IsolatedAtomPotential(const std::map<Species, double> &isolated_atom_energies, bool error_on_unknown) {
        isolated_energies_ = isolated_atom_energies;
        error_on_unknown_species_ = error_on_unknown;
    }

    DataNode IsolatedAtomPotential::serialize() {
        return DataNode::object(){
            {"error_on_unknown", error_on_unknown_species_},
            {"energies", isolated_energies_}
        };
    }

    Predictions IsolatedAtomPotential::predict(const AtomicStructure& structure) {

        double result = 0;

        for (const auto &atom: structure) {
            if (isolated_energies_.contains(atom.species())) {
                result += isolated_energies_[atom.species()];
            } else {
                if (error_on_unknown_species_) {
                    JGAP_LOG -> error("Unknown isolated_atom energy for " + atom.species(),true);
                }
            }
        }

        return Predictions{.energy = result};
    }

    void IsolatedAtomPotential::tabulate(TabulationData &table) {
        for (const auto &[species, energy]: isolated_energies_) {
            if (table.isolatedEnergies.contains(species)) {
                JGAP_LOG_WARN("Conflicting isolated atom energies for {}", species);
            } else {
                table.isolatedEnergies[species] = 0;
            }
            table.isolatedEnergies[species] += energy;
        }
    }
}
