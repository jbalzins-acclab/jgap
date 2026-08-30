#include "PerConfigTypeRegularizationRules.hpp"
#include <sstream>
#include "../../../io/log/CurrentLogger.hpp"
#include "jgap/core/atomic/Atoms.hpp"

namespace jgap {
    PerConfigTypeRegularizationRules::PerConfigTypeRegularizationRules(
        PerConfigTypeSigmas default_sigmas,
        std::map<std::string, PerConfigTypeSigmas> exact_config_type_sigmas,
        std::map<std::string, PerConfigTypeSigmas> config_type_contains_sigmas
    ) :
        defaults(std::move(default_sigmas)),
        exact_config_type_sigmas(std::move(exact_config_type_sigmas)),
        config_type_contains_sigmas(std::move(config_type_contains_sigmas)) {}

    PerConfigTypeRegularizationRules::PerConfigTypeRegularizationRules(
        PerConfigTypeSigmas default_sigmas,
        const std::string& config_string
    ) :
        defaults(std::move(default_sigmas)) {

        std::stringstream ss(config_string);
        std::string segment;
        std::vector<std::string> segments;
        while (std::getline(ss, segment, ':')) {
            segments.push_back(segment);
        }

        if (segments.size() % 5 != 0) {
            JGAP_LOG_AND_THROW("Invalid config string format: {}", config_string);
        }

        for (size_t i = 0; i < segments.size(); i += 5) {
            std::string config_type = segments[i];
            try {
                Real e = std::stod(segments[i + 1]);
                Real f = std::stod(segments[i + 2]);
                Real v = std::stod(segments[i + 3]);
                Real h = std::stod(segments[i + 4]);
                if (h != 0.0) {
                    JGAP_LOG_WARN(
                        "Hessian regularization is not supported, but a non-zero value was provided for config_type {}",
                        config_type
                    );
                }
                exact_config_type_sigmas.insert({config_type, PerConfigTypeSigmas(e, f, v)});
            } catch (const std::invalid_argument& e) {
                JGAP_LOG_AND_THROW("Invalid number in config string: {}", config_string);
            }
        }
    }

    Regularization PerConfigTypeRegularizationRules::determine(const Atoms& atoms) const {
        const std::string ct = atoms.getConfigType().value_or("");

        const PerConfigTypeSigmas* matched_sigmas = nullptr;

        if (auto it = exact_config_type_sigmas.find(ct); it != exact_config_type_sigmas.end()) {
            matched_sigmas = &it->second;
        } else {
            std::string longest_match_key;
            for (const auto& [key, val]: config_type_contains_sigmas) {
                if (ct.find(key) != std::string::npos) {
                    if (key.length() > longest_match_key.length()) {
                        longest_match_key = key;
                    }
                }
            }
            if (!longest_match_key.empty()) {
                matched_sigmas = &config_type_contains_sigmas.at(longest_match_key);
            }
        }

        const auto& s = (matched_sigmas != nullptr) ? *matched_sigmas : defaults;

        Regularization sigmas;
        sigmas.energy = s.energy;
        sigmas.forces = std::vector<Vector3>(atoms.nAtoms(), s.force);
        sigmas.virials = s.virials;
        return sigmas;
    }
}