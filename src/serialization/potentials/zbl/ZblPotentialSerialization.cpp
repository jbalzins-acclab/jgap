#include "ZblPotentialSerialization.hpp"
#include "core/potentials/zbl/ZblPotential.hpp"
#include "serialization/SerializationNode.hpp"
#include "io/log/CurrentLogger.hpp"

namespace jgap {

    bool ZblPotentialSerialization::serialize(const ValuePtr<Potential>& obj, SerializationNode& node) const {
        if (auto derived = obj.as<ZblPotential>()) {
            node.writeAttribute("name", "ZblPotential");
            node.writeAttribute("cutoff", derived->getCutoff());
            node.writeAttribute("cutoff_transition_width", derived->getCutoffTransitionWidth());

            // z1_z2 and a_inverse are deduced from the species, so only the coefficients are stored.
            auto coefficients_group = node.createGroup("coefficients");
            int i = 0;
            for (const auto& [species_set, coeffs] : derived->getCoefficients()) {
                auto pair_group = coefficients_group.createGroup(std::to_string(i++));

                std::vector<std::string> species_symbols;
                for (const auto& s : species_set.getNodes()) {
                    species_symbols.push_back(s.symbol());
                }
                pair_group.writeAttribute("species_set", species_symbols);
                pair_group.writeDataSet("coeffs", coeffs);
            }
            return true;
        }
        return false;
    }

    ValuePtr<Potential> ZblPotentialSerialization::deserialize(const SerializationNode& node) const {
        if (node.readOptionalAttribute<std::string>("name") != "ZblPotential") {
            return nullptr;
        }

        auto cutoff = node.readAttribute<Real>("cutoff");
        auto cutoff_transition_width = node.readAttribute<Real>("cutoff_transition_width");

        auto coefficients_group_opt = node.getGroup("coefficients");
        if (!coefficients_group_opt) {
            JGAP_LOG_AND_THROW("Missing 'coefficients' group in ZblPotential serialization");
        }
        const auto& coefficients_group = coefficients_group_opt.value();

        std::map<SpeciesSet<2, Symmetric>, std::array<Real, 6>> coefficients;
        for (const auto& group_name : coefficients_group.getChildNames()) {
            auto pair_group_opt = coefficients_group.getGroup(group_name);
            if (!pair_group_opt) {
                JGAP_LOG_AND_THROW("Missing coefficient group in ZblPotential serialization");
            }
            const auto& pair_group = pair_group_opt.value();

            auto species_symbols = pair_group.readAttribute<std::vector<std::string>>("species_set");
            if (species_symbols.size() != 2) {
                JGAP_LOG_AND_THROW("Expected 2 species in species_set for ZblPotential");
            }
            SpeciesSet<2, Symmetric> pair{Species(species_symbols[0]), Species(species_symbols[1])};

            coefficients[pair] = pair_group.readDataSet<std::array<Real, 6>>("coeffs");
        }

        return ValuePtr<Potential>(ZblPotential(coefficients, cutoff, cutoff_transition_width));
    }

    REGISTER_SERIALIZATION(ZblPotentialSerialization, Potential);
}
