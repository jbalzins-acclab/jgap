#include "core/matrices/regularization/PerConfigTypeRegularizationRules.hpp"

#include "io/log/CurrentLogger.hpp"
#include "utils/Utils.hpp"

namespace jgap {
    std::shared_ptr<PerConfigTypeRegularizationRules> PerConfigTypeRegularizationRules::fromDataNode(
                                                        const DataNode &params) {

        double e = require(params, Conventions::ENERGY_PER_ATOM);
        double f = params.getOrDefault(Conventions::FORCE_COMPONENT, e * 50.0);

        double vIso, vAniso;
        if (params.contains(Conventions::VIRIAL_ISO_PER_ATOM)
            || params.contains(Conventions::VIRIAL_ANISO_PER_ATOM)) {

            vIso = require(params, Conventions::VIRIAL_ISO_PER_ATOM);
            vAniso = require(params, Conventions::VIRIAL_ANISO_PER_ATOM);

        } else if (params.contains(Conventions::VIRIAL_PER_ATOM)) {
            vIso = vAniso = require(params, Conventions::VIRIAL_PER_ATOM);
        } else {
            vIso = vAniso = e * 100;
        }

        std::map<std::string, double> exact;
        if (params.contains("exact")) {
            for (const auto &[ct, multiplier]: require(params, "exact").asObject()) {
                if (exact.contains(ct)) {
                    JGAP_LOG_AND_THROW("Multiple 'exact' rules for config_type={}", ct);
                }
                exact[ct] = multiplier;
            }
        }
        std::map<std::string, double> contains;
        if (params.contains("contains")) {
            for (const auto &[ct, multiplier]: require(params, "contains").asObject()) {
                if (contains.contains(ct)) {
                    JGAP_LOG_AND_THROW("Multiple 'contains' rules for config_type={}", ct);
                }
                contains[ct] = multiplier;
            }
        }

        return std::make_shared<PerConfigTypeRegularizationRules>(e, f, vIso, vAniso, exact, contains);
    }

    PerConfigTypeRegularizationRules::PerConfigTypeRegularizationRules(
        double energyPerAtom, double forceComponent, double virialIsoPerAtom, double virialAnisoPerAtom,
        const std::map<std::string, double> &exact, const std::map<std::string, double> &contains)
        : _defaultEnergyPerAtom(energyPerAtom), _defaultForceComponent(forceComponent),
          _defaultVirialIsoPerAtom(virialIsoPerAtom), _defaultVirialAnisoPerAtom(virialAnisoPerAtom),
          _exact(exact), _contains(contains) {
    }

    void PerConfigTypeRegularizationRules::fillSigmas(std::vector<AtomicStructure> &structures) {
        for (auto& structure: structures) {
            fillSigmas(structure);
        }
    }

    void PerConfigTypeRegularizationRules::fillSigmas(AtomicStructure &structure) {
        double multiplier = 1.0;

        if (structure.properties.contains("config_type")) {
            const auto &ct = structure.properties.at("config_type");

            if (_exact.contains(ct)) {
                multiplier = _exact.at(ct);
            } else {
                for (const auto &[subString, multiplierPerStr]: _contains) {
                    if (ct.contains(subString)) {
                        multiplier = multiplierPerStr;
                        break;
                    }
                }
            }
        }

        RegularizationRules::fillSigmas(
            structure, multiplier * _defaultEnergyPerAtom,
            Vector3{_defaultForceComponent, _defaultForceComponent, _defaultForceComponent} * multiplier,
            multiplier * _defaultVirialIsoPerAtom, multiplier * _defaultVirialAnisoPerAtom
            );
    }
}
