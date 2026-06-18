#include "IsolatedAtomPotentialSerialization.hpp"
#include "core/potentials/isolated/IsolatedAtomPotential.hpp"
#include "serialization/SerializationNode.hpp"
#include "io/log/CurrentLogger.hpp"

namespace jgap {

    bool IsolatedAtomPotentialSerialization::serialize(const ValuePtr<Potential>& obj, SerializationNode& node) const {
        if (auto derived = obj.as<IsolatedAtomPotential>()) {
            node.writeAttribute("name", "IsolatedAtomPotential");

            auto energies_group = node.createGroup("isolated_energies");
            for (const auto& [species, energy] : derived->getIsolatedEnergies()) {
                energies_group.writeAttribute(species.symbol(), energy);
            }
            return true;
        }
        return false;
    }

    ValuePtr<Potential> IsolatedAtomPotentialSerialization::deserialize(const SerializationNode& node) const {
        if (node.readOptionalAttribute<std::string>("name") != "IsolatedAtomPotential") {
            return nullptr;
        }

        std::map<Species, Real> isolated_energies;
        auto energies_group_opt = node.getGroup("isolated_energies");
        if (!energies_group_opt) {
            JGAP_LOG_AND_THROW("Missing 'isolated_energies' group in IsolatedAtomPotential serialization");
        }

        const auto& energies_group = energies_group_opt.value();
        for (const auto& symbol : energies_group.getAttributeNames()) {
            auto energy = energies_group.readAttribute<Real>(symbol);
            isolated_energies.emplace(Species(symbol), energy);
        }

        return ValuePtr<Potential>(IsolatedAtomPotential(isolated_energies));
    }

    REGISTER_SERIALIZATION(IsolatedAtomPotentialSerialization, Potential);
}