#include "ScreenedCoulombPotentialSerialization.hpp"
#include "jgap/core/potentials/coulomb/ScreenedCoulombPotential.hpp"
#include "jgap/core/io/log/CurrentLogger.hpp"
#include "jgap/serialization/SerializationNode.hpp"

#include "jgap/core/atomic/species/composition/Species2Sorted.hpp"

namespace jgap {

    bool ScreenedCoulombPotentialSerialization::serialize(const ValuePtr<Potential>& obj, SerializationNode& node) const {
        if (auto derived = obj.as<ScreenedCoulombPotential>()) {
            node.writeAttribute("name", "ScreenedCoulombPotential");
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

    ValuePtr<Potential> ScreenedCoulombPotentialSerialization::deserialize(const SerializationNode& node) const {
        if (node.readOptionalStringAttribute("name") != "ScreenedCoulombPotential") {
            return nullptr;
        }

        auto cutoff = node.readDoubleAttribute("cutoff");
        auto cutoff_transition_width = node.readDoubleAttribute("cutoff_transition_width");

        auto coefficients_group_opt = node.getGroup("coefficients");
        if (!coefficients_group_opt) {
            JGAP_LOG_AND_THROW("Missing 'coefficients' group in ScreenedCoulombPotential serialization");
        }
        const auto& coefficients_group = coefficients_group_opt.value();

        std::map<Species2Sorted, std::array<Real, 6>> coefficients;
        for (const auto& group_name: coefficients_group.getChildNames()) {
            auto pair_group_opt = coefficients_group.getGroup(group_name);
            if (!pair_group_opt) {
                JGAP_LOG_AND_THROW("Missing coefficient group in ScreenedCoulombPotential serialization");
            }
            const auto& pair_group = pair_group_opt.value();

            auto species_encoded = pair_group.readStringAttribute("species_set");
            Species2Sorted pair(species_encoded);

            coefficients[pair] = pair_group.readDoubleArrayDataSet<6>("coeffs");
        }

        return ValuePtr<Potential>(ScreenedCoulombPotential(coefficients, cutoff, cutoff_transition_width));
    }

    REGISTER_SERIALIZATION(ScreenedCoulombPotentialSerialization, Potential);
}
