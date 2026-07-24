#include "ZblPotentialSerialization.hpp"
#include "jgap/core/potentials/zbl/ZblPotential.hpp"
#include "jgap/io/log/CurrentLogger.hpp"
#include "jgap/serialization/SerializationNode.hpp"

#include "jgap/core/atomic/species/composition/Species2Sorted.hpp"

namespace jgap {

    bool ZblPotentialSerialization::serialize(const ValuePtr<Potential>& obj, SerializationNode& node) const {
        if (auto derived = obj.as<ZblPotential>()) {
            node.writeAttribute("name", "ZblPotential");
            node.writeAttribute("cutoff", derived->getCutoff());
            node.writeAttribute("cutoff_transition_width", derived->getCutoffTransitionWidth());

            // z1_z2 and a_inverse are deduced from the species, so only the coefficients are stored.
            auto coefficients_group = node.createGroup("coefficients");
            int i = 0;
            for (const auto& [species_set, coeffs]: derived->getCoefficients()) {
                auto pair_group = coefficients_group.createGroup(std::to_string(i++));

                pair_group.writeAttribute("species_set", species_set.toString());
                pair_group.writeDataSet("coeffs", coeffs);
            }
            return true;
        }
        return false;
    }

    ValuePtr<Potential> ZblPotentialSerialization::deserialize(const SerializationNode& node) const {
        if (node.readOptionalStringAttribute("name") != "ZblPotential") {
            return nullptr;
        }

        auto cutoff = node.readDoubleAttribute("cutoff");
        auto cutoff_transition_width = node.readDoubleAttribute("cutoff_transition_width");

        auto coefficients_group_opt = node.getGroup("coefficients");
        if (!coefficients_group_opt) {
            JGAP_LOG_AND_THROW("Missing 'coefficients' group in ZblPotential serialization");
        }
        const auto& coefficients_group = coefficients_group_opt.value();

        std::map<Species2Sorted, std::array<Real, 6>> coefficients;
        for (const auto& group_name: coefficients_group.getChildNames()) {
            auto pair_group_opt = coefficients_group.getGroup(group_name);
            if (!pair_group_opt) {
                JGAP_LOG_AND_THROW("Missing coefficient group in ZblPotential serialization");
            }
            const auto& pair_group = pair_group_opt.value();

            auto species_encoded = pair_group.readStringAttribute("species_set");
            Species2Sorted pair(species_encoded);

            coefficients[pair] = pair_group.readDoubleArrayDataSet<6>("coeffs");
        }

        return ValuePtr<Potential>(ZblPotential(coefficients, cutoff, cutoff_transition_width));
    }

    REGISTER_SERIALIZATION(ZblPotentialSerialization, Potential);
}
