//
// Created by Jegors Balzins on 8.8.2025.
//

#include "../../../../include/core/matrices/sigmas/PerConfigTypeRegularizationRules.hpp"

namespace jgap {
    PerConfigTypeRegularizationRules::PerConfigTypeRegularizationRules(const nlohmann::json &params) {
        _E = params["E"];
        _V = params["V"];
        _F = params["F"];

        // WARN: order sensitive!
        _multipliersPerKeyWord = {};
        for (const auto &perKeyWord : params["per_keyword"]) {
            _multipliersPerKeyWord.push_back({
                perKeyWord["contains"].get<string>(),
                perKeyWord["multiplier"].get<double>(),
            });
        }
        /*
         {
            "type": "per_ct",
            "E": 0.001,
            "F": 0.05,
            "V": 0.1,
            "per_keyword": [
                {
                    "contains": "isolated_atom",
                    "multiplier": 0.001
                },
                {
                    "contains": "liquid_high",
                    "multiplier": 10.0
                },
                {
                    "contains": "liquid",
                    "multiplier": 5.0
                },
                {
                    "contains": "short",
                    "multiplier": 5.0
                }
            ]
         }
         */
    }

    void PerConfigTypeRegularizationRules::fillSigmas(AtomicStructure &structure) {
        const double mul = structure.configType.transform([&](string ct) -> double {
            for (const auto &[keyWord, multiplier] : _multipliersPerKeyWord) {
                if (ct.contains(keyWord)) {
                    // structure.virials.reset();
                    return multiplier;
                }
            }
            return 1.0;
        }).value_or(1.0);

        const double E = _E * mul;
        const double F = _F * mul;
        const double V = _V * mul;

        if (!structure.energySigmaInverse.has_value()) {
            structure.energySigmaInverse = 1.0 / (E * pow(structure.size(), 0.5));
        }

        const double dF = 1.0 / F;
        if (!structure.forceSigmasInverse.has_value()) {
            structure.forceSigmasInverse = vector(structure.size(),  Vector3{dF, dF, dF});
        }

        const double dV = 1.0 / (V * pow(structure.size(), 0.5));
        if (!structure.virialSigmasInverse.has_value()) {
            structure.virialSigmasInverse = {
                Vector3{dV, dV, dV},
                Vector3{dV, dV, dV},
                Vector3{dV, dV, dV}
            };
        }
    }
}
