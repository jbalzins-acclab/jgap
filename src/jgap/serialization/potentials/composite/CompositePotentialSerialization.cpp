#include "CompositePotentialSerialization.hpp"
#include "jgap/core/potentials/CompositePotential.hpp"
#include "jgap/serialization/SerializationNode.hpp"
#include "../../../core/io/log/CurrentLogger.hpp"

namespace jgap {

    bool CompositePotentialSerialization::serialize(const ValuePtr<Potential>& obj, SerializationNode& node) const {
        const auto derived = obj.as<CompositePotential>();
        if (!derived) {
            return false;
        }

        node.writeAttribute("name", "CompositePotential");

        // Keys are arbitrary strings, so they are stored as an attribute on numbered subgroups
        // rather than used as (potentially invalid) HDF5 group names.
        auto potentials_group = node.createGroup("potentials");
        int i = 0;
        for (const auto& [key, potential] : derived->getPotentials()) {
            auto entry_group = potentials_group.createGroup(std::to_string(i++));
            entry_group.writeAttribute("key", key);
            auto potential_group = entry_group.createGroup("potential");
            SerializationRegistry<Potential>::serialize(potential, potential_group);
        }

        return true;
    }

    ValuePtr<Potential> CompositePotentialSerialization::deserialize(const SerializationNode& node) const {
        if (node.readOptionalStringAttribute("name") != "CompositePotential") {
            return nullptr;
        }

        auto potentials_group_opt = node.getGroup("potentials");
        if (!potentials_group_opt) {
            JGAP_LOG_AND_THROW("Missing 'potentials' group in CompositePotential serialization");
        }
        const auto& potentials_group = potentials_group_opt.value();

        std::map<std::string, ValuePtr<Potential>> potentials;
        for (const auto& group_name : potentials_group.getChildNames()) {
            auto entry_group_opt = potentials_group.getGroup(group_name);
            if (!entry_group_opt) {
                JGAP_LOG_AND_THROW("Missing entry group in CompositePotential serialization");
            }
            const auto& entry_group = entry_group_opt.value();

            auto key = entry_group.readStringAttribute("key");
            auto potential_group_opt = entry_group.getGroup("potential");
            if (!potential_group_opt) {
                JGAP_LOG_AND_THROW("Missing 'potential' group in CompositePotential serialization");
            }
            potentials[key] = SerializationRegistry<Potential>::deserialize(potential_group_opt.value());
        }

        return ValuePtr<Potential>(CompositePotential(std::move(potentials)));
    }

    REGISTER_SERIALIZATION(CompositePotentialSerialization, Potential);
}
