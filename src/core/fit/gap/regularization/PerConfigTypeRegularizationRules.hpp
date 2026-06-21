#ifndef JGAP_PERCONFIGTYPEREGULARIZATIONRULES_HPP
#define JGAP_PERCONFIGTYPEREGULARIZATIONRULES_HPP

#include <optional>
#include "RegularizationRules.hpp"
#include "ConfigSigmas.hpp"
#include "../../../../serialization/SerializationRegistry.hpp"

namespace jgap {
    class PerConfigTypeRegularizationRules : public RegularizationRules {
    public:
        PerConfigTypeRegularizationRules(Real energy_sigma_per_atom,
                                         Real force_component_sigma,
                                         Real virials_iso_sigma_per_atom,
                                         Real virials_aniso_sigmas_per_atom,
                                         const std::map<std::string, Real> &exact,
                                         const std::map<std::string, Real> &contains
        );

        PerConfigTypeRegularizationRules(const ConfigSigmas &default_sigmas,
                                         const std::map<std::string, ConfigSigmas>& exact_config_type_sigmas,
                                         const std::map<std::string, ConfigSigmas>& contains_config_type_sigmas);

        PerConfigTypeRegularizationRules(ConfigSigmas default_sigmas,
                                         const std::string& config_string);

        void fillSigmas(Regularization &sigmas, const Atoms &atoms) const override;

        PerConfigTypeRegularizationRules* clone() const override {
            return new PerConfigTypeRegularizationRules(*this);
        }

    private:
        ConfigSigmas defaults;

        std::map<std::string, Real> exact_multiplier;
        std::map<std::string, Real> contains_multiplier; // order-sensitive
        std::map<std::string, ConfigSigmas> exact_config_type_sigmas;
        std::map<std::string, ConfigSigmas> contains_config_type_sigmas;
    };
}

#endif