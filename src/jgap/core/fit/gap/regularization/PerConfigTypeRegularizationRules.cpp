#include "PerConfigTypeRegularizationRules.hpp"
#include "jgap/core/atomic/Atoms.hpp"
#include <cmath>
#include <sstream>
#include <vector>
#include "jgap/io/log/CurrentLogger.hpp"

namespace jgap {
    PerConfigTypeRegularizationRules::PerConfigTypeRegularizationRules(Real energy_sigma_per_atom,
                                                                       Real force_component_sigma,
                                                                       Real virials_iso_sigma_per_atom,
                                                                       Real virials_aniso_sigmas_per_atom,
                                                                       const std::map<std::string, Real> &exact,
                                                                       const std::map<std::string, Real> &contains)
        : defaults(energy_sigma_per_atom, force_component_sigma, virials_iso_sigma_per_atom, virials_aniso_sigmas_per_atom),
          exact_multiplier(exact),
          contains_multiplier(contains) {
    }

    PerConfigTypeRegularizationRules::PerConfigTypeRegularizationRules(const ConfigSigmas &default_sigmas,
        const std::map<std::string, ConfigSigmas> &exact_config_type_sigmas,
        const std::map<std::string, ConfigSigmas> &contains_config_type_sigmas)
        : defaults(default_sigmas),
          exact_config_type_sigmas(exact_config_type_sigmas),
          contains_config_type_sigmas(contains_config_type_sigmas) {
    }

    PerConfigTypeRegularizationRules::PerConfigTypeRegularizationRules(ConfigSigmas default_sigmas,
                                                                       const std::string& config_string)
        : defaults(default_sigmas) {

        std::stringstream ss(config_string);
        std::string segment;
        std::vector<std::string> segments;
        while(std::getline(ss, segment, ':')) {
            segments.push_back(segment);
        }

        if (segments.size() % 5 != 0) {
            JGAP_LOG_AND_THROW("Invalid config string format: {}", config_string);
        }

        for (size_t i = 0; i < segments.size(); i += 5) {
            std::string config_type = segments[i];
            try {
                Real e = std::stod(segments[i+1]);
                Real f = std::stod(segments[i+2]);
                Real v = std::stod(segments[i+3]);
                Real h = std::stod(segments[i+4]);
                if (h != 0.0) {
                    JGAP_LOG_WARN("Hessian regularization is not supported, "
                                  "but a non-zero value was provided for config_type {}", config_type);
                }
                exact_config_type_sigmas.insert({config_type, ConfigSigmas(e, f, v)});
            } catch (const std::invalid_argument& e) {
                JGAP_LOG_AND_THROW("Invalid number in config string: {}", config_string);
            }
        }
    }

    void PerConfigTypeRegularizationRules::fillSigmas(Regularization& sigmas, const Atoms& atoms) const {
        const std::string ct = atoms.getConfigType().value_or("");

        if (exact_config_type_sigmas.contains(ct)) {
            const auto& s = exact_config_type_sigmas.at(ct);
            sigmas.energy = s.energy;
            sigmas.forces = std::vector<Vector3>(atoms.nAtoms(), s.force);
            sigmas.virials = s.virials;
            return;
        }

        if (!contains_config_type_sigmas.empty()) {
            std::string longest_match_key;
            for (const auto& [key, value] : contains_config_type_sigmas) {
                if (ct.find(key) != std::string::npos) {
                    if (key.length() > longest_match_key.length()) {
                        longest_match_key = key;
                    }
                }
            }
            if (!longest_match_key.empty()) {
                const auto& s = contains_config_type_sigmas.at(longest_match_key);
                sigmas.energy = s.energy;
                sigmas.forces = std::vector<Vector3>(atoms.nAtoms(), s.force);
                sigmas.virials = s.virials;
                return;
            }
        }

        Real multiplier = 1.0;
        if (exact_multiplier.contains(ct)) {
            multiplier = exact_multiplier.at(ct);
        } else {
            std::string longest_match_key;
            for (const auto& [key, value] : contains_multiplier) {
                if (ct.find(key) != std::string::npos) {
                    if (key.length() > longest_match_key.length()) {
                        longest_match_key = key;
                    }
                }
            }

            if (!longest_match_key.empty()) {
                multiplier = contains_multiplier.at(longest_match_key);
            }
        }

        sigmas.energy = defaults.energy * multiplier * pow(atoms.nAtoms(), 0.5);
        sigmas.virials = defaults.virials * multiplier * pow(atoms.nAtoms(), 0.5);

        sigmas.forces = std::vector<Vector3>(atoms.nAtoms());
        for (int i = 0; i < atoms.nAtoms(); i++) {
            (*sigmas.forces)[i] = defaults.force * multiplier;
        }
    }
}