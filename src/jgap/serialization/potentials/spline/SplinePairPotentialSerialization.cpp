#include "SplinePairPotentialSerialization.hpp"
#include "jgap/core/potentials/spline/SplinePairPotential.hpp"
#include "jgap/io/log/CurrentLogger.hpp"
#include "jgap/serialization/SerializationNode.hpp"

namespace jgap {

    bool SplinePairPotentialSerialization::serialize(const ValuePtr<Potential>& obj, SerializationNode& node) const {
        if (auto derived = obj.as<SplinePairPotential>()) {
            node.writeAttribute("name", "SplinePairPotential");

            auto interpolators_group = node.createGroup("interpolators");
            int i = 0;
            for (const auto& [species_set, spline]: derived->getInterpolators()) {
                auto spline_group = interpolators_group.createGroup(std::to_string(i++));

                spline_group.writeAttribute("species_set", species_set.toString());

                spline_group.writeDataSet("r_vec", spline.getRVec());
                spline_group.writeDataSet("energies", spline.getEnergies());
            }
            return true;
        }
        return false;
    }

    ValuePtr<Potential> SplinePairPotentialSerialization::deserialize(const SerializationNode& node) const {
        if (node.readOptionalStringAttribute("name") != "SplinePairPotential") {
            return nullptr;
        }

        auto potential = std::make_unique<SplinePairPotential>();

        auto interpolators_group_opt = node.getGroup("interpolators");
        if (!interpolators_group_opt) {
            JGAP_LOG_AND_THROW("Missing 'interpolators' group in SplinePairPotential serialization");
        }
        const auto& interpolators_group = interpolators_group_opt.value();

        for (const auto& group_name: interpolators_group.getChildNames()) {
            auto spline_group_opt = interpolators_group.getGroup(group_name);
            if (!spline_group_opt) {
                JGAP_LOG_AND_THROW("Missing interpolator group in SplinePairPotential serialization");
            }
            const auto& spline_group = spline_group_opt.value();

            auto species_encoded = spline_group.readStringAttribute("species_set");
            Species2Sorted species_set(species_encoded);
            Species s1 = species_set.nodes[0];
            Species s2 = species_set.nodes[1];

            auto r_vec = spline_group.readRealVectorDataSet("r_vec");
            auto energies = spline_group.readRealVectorDataSet("energies");

            potential->extend(s1, s2, r_vec, energies);
        }

        return ValuePtr<Potential>(std::move(potential));
    }

    REGISTER_SERIALIZATION(SplinePairPotentialSerialization, Potential);
}
