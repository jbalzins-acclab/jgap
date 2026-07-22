#include "GapPotentialSerialization.hpp"
#include "core/potentials/gap/GapPotential.hpp"
#include "core/potentials/gap/component/GapComponent.hpp"
#include "io/log/CurrentLogger.hpp"
#include "serialization/SerializationNode.hpp"
#include "serialization/SerializationRegistry.hpp"

#include <algorithm>
#include <string>

namespace jgap {

    bool GapPotentialSerialization::serialize(const ValuePtr<Potential>& obj, SerializationNode& node) const {
        const auto derived = obj.as<GapPotential>();
        if (!derived) {
            return false;
        }

        node.writeAttribute("name", "GapPotential");

        if (derived->optional_external_potential.get() != nullptr) {
            auto external_group = node.createGroup("external_potential");
            SerializationRegistry<Potential>::serialize(derived->optional_external_potential, external_group);
        }

        auto components_group = node.createGroup("components");
        int i = 0;
        for (const auto& component: derived->getComponents()) {
            auto component_group = components_group.createGroup(std::to_string(i++));
            SerializationRegistry<GapComponent>::serialize(component, component_group);
        }

        return true;
    }

    ValuePtr<Potential> GapPotentialSerialization::deserialize(const SerializationNode& node) const {
        if (node.readOptionalAttribute<std::string>("name") != "GapPotential") {
            return nullptr;
        }

        auto potential = std::make_unique<GapPotential>();

        if (auto external_group_opt = node.getGroup("external_potential")) {
            potential->optional_external_potential =
                SerializationRegistry<Potential>::deserialize(external_group_opt.value());
        }

        auto components_group_opt = node.getGroup("components");
        if (!components_group_opt) {
            JGAP_LOG_AND_THROW("Missing 'components' group in GapPotential serialization");
        }
        const auto& components_group = components_group_opt.value();

        // Child group names are numeric indices; sort numerically to preserve component order.
        auto child_names = components_group.getChildNames();
        std::sort(child_names.begin(), child_names.end(),
                  [](const std::string& a, const std::string& b) { return std::stoi(a) < std::stoi(b); });

        for (const auto& group_name: child_names) {
            auto component_group_opt = components_group.getGroup(group_name);
            if (!component_group_opt) {
                JGAP_LOG_AND_THROW("Missing component group in GapPotential serialization");
            }
            potential->addComponent(SerializationRegistry<GapComponent>::deserialize(component_group_opt.value()));
        }

        return ValuePtr<Potential>(std::move(potential));
    }

    REGISTER_SERIALIZATION(GapPotentialSerialization, Potential);
}
