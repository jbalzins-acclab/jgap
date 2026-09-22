#ifndef JGAP_PERCONFIGTYPEREGULARIZATIONRULES_HPP
#define JGAP_PERCONFIGTYPEREGULARIZATIONRULES_HPP

#include <map>
#include <string>
#include "PerConfigTypeSigmas.hpp"
#include "RegularizationRules.hpp"

namespace jgap {
    class PerConfigTypeRegularizationRules : public RegularizationRules {
    public:
        explicit PerConfigTypeRegularizationRules(
            PerConfigTypeSigmas default_sigmas,
            std::map<std::string, PerConfigTypeSigmas> exact_config_type_sigmas = {},
            std::map<std::string, PerConfigTypeSigmas> config_type_contains_sigmas = {}
        );

        PerConfigTypeRegularizationRules(PerConfigTypeSigmas default_sigmas, const std::string& config_string);

        Regularization determine(const Atoms& atoms) const override;

        PerConfigTypeRegularizationRules* clone() const override { return new PerConfigTypeRegularizationRules(*this); }

        const PerConfigTypeSigmas& getDefaults() const { return defaults; }
        const std::map<std::string, PerConfigTypeSigmas>& getExactConfigTypeSigmas() const {
            return exact_config_type_sigmas;
        }
        const std::map<std::string, PerConfigTypeSigmas>& getConfigTypeContainsSigmas() const {
            return config_type_contains_sigmas;
        }

    private:
        PerConfigTypeSigmas defaults;
        std::map<std::string, PerConfigTypeSigmas> exact_config_type_sigmas;
        std::map<std::string, PerConfigTypeSigmas> config_type_contains_sigmas;
    };
}

#endif
