#include "SplinePairPotentialSerialization.hpp"
#include "core/potentials/spline/SplinePairPotential.hpp"
#include "serialization/SerializationNode.hpp"
#include "io/log/CurrentLogger.hpp"

namespace jgap {

    bool SplinePairPotentialSerialization::serialize(const ValuePtr<Potential>& obj, SerializationNode& node) const {
        if (auto derived = obj.as<SplinePairPotential>()) {
            node.writeAttribute("name", "SplinePairPotential");

            auto interpolators_group = node.createGroup("interpolators");
            int i = 0;
            for (const auto& [species_set, spline] : derived->getInterpolators()) {
                auto spline_group = interpolators_group.createGroup(std::to_string(i++));

                std::vector<std::string> species_symbols;
                for (const auto& s : species_set.getNodes()) {
                    species_symbols.push_back(s.symbol());
                }
                spline_group.writeAttribute("species_set", species_symbols);

                spline_group.writeDataSet("r_vec", spline.getRVec());
                spline_group.writeDataSet("energies", spline.getEnergies());
            }
            return true;
        }
        return false;
    }

    ValuePtr<Potential> SplinePairPotentialSerialization::deserialize(const SerializationNode& node) const {
        if (node.readOptionalAttribute<std::string>("name") != "SplinePairPotential") {
            return nullptr;
        }

        auto potential = std::make_unique<SplinePairPotential>();

        auto interpolators_group_opt = node.getGroup("interpolators");
        if (!interpolators_group_opt) {
            JGAP_LOG_AND_THROW("Missing 'interpolators' group in SplinePairPotential serialization");
        }
        const auto& interpolators_group = interpolators_group_opt.value();

        for (const auto& group_name : interpolators_group.getChildNames()) {
            auto spline_group_opt = interpolators_group.getGroup(group_name);
            if (!spline_group_opt) {
                JGAP_LOG_AND_THROW("Missing interpolator group in SplinePairPotential serialization");
            }
            const auto& spline_group = spline_group_opt.value();

            auto species_symbols = spline_group.readAttribute<std::vector<std::string>>("species_set");
            if (species_symbols.size() != 2) {
                 JGAP_LOG_AND_THROW("Expected 2 species in species_set for SplinePairPotential");
            }
            Species s1(species_symbols[0]);
            Species s2(species_symbols[1]);

            auto r_vec = spline_group.readDataSet<std::vector<Real>>("r_vec");
            auto energies = spline_group.readDataSet<std::vector<Real>>("energies");

            potential->extend(s1, s2, r_vec, energies);
        }

        return ValuePtr<Potential>(std::move(potential));
    }

    REGISTER_SERIALIZATION(SplinePairPotentialSerialization, Potential);
}