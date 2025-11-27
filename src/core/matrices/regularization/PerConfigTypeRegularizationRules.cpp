#include "core/matrices/regularization/PerConfigTypeRegularizationRules.hpp"

#include "utils/Utils.hpp"

namespace jgap {
    PerConfigTypeRegularizationRules::PerConfigTypeRegularizationRules(const nlohmann::json &params) {
        _energy = require(params, "E");
        _force = require(params, "F");

        // Virial split handling with backward compatibility
        if (params.contains("V_iso") || params.contains("V_aniso")) {
            const double v_iso = params.contains("V_iso") ? params["V_iso"].get<double>() :
                                 (params.contains("V") ? params["V"].get<double>() : 0.0);
            const double v_aniso = params.contains("V_aniso") ? params["V_aniso"].get<double>() : v_iso;
            _virialIso = v_iso;
            _virialAniso = v_aniso;
        } else {
            const double v = params.contains("V") ? params["V"].get<double>() : 0.0;
            _virialIso = v;
            _virialAniso = v;
        }

        // WARN: order sensitive within the same rule type (exact/contains)!
        _rules.clear();
        if (params.contains("per_keyword")) {
            for (const auto &perKeyWord : params["per_keyword"]) {
                Rule rule{};
                rule.multiplier = perKeyWord["multiplier"];
                if (perKeyWord.contains("is")) {
                    rule.exact = true;
                    rule.key = perKeyWord["is"];
                } else {
                    rule.exact = false;
                    rule.key = perKeyWord["contains"];
                }
                _rules.push_back(rule);
            }
        }
        /*
         {
            "type": "per_config_type",
            "E": 0.001,
            "F": 0.05,
            // Either legacy single virial value V or split values V_iso/V_aniso
            "V": 0.1,
            // or
            // "V_iso": 0.1,
            // "V_aniso": 0.2,
            "per_keyword": [
                { "is": "isolated_atom", "multiplier": 0.001 },
                { "contains": "liquid_high", "multiplier": 10.0 },
                { "contains": "liquid", "multiplier": 5.0 },
                { "contains": "short", "multiplier": 5.0 }
            ]
         }
         */
    }

    void PerConfigTypeRegularizationRules::fillSigmas(AtomicStructure &structure) {
        double mul = 1.0;
        if (structure.properties.contains("config_type")) {
            const auto &ct = structure.properties.at("config_type");

            // Two-pass priority: track last matching exact rule; otherwise last matching contains rule
            optional<double> mulExact;
            optional<double> mulContains;
            for (const auto &rule : _rules) {
                if (rule.exact) {
                    if (ct == rule.key) mulExact = rule.multiplier; // order sensitive: last wins
                } else {
                    if (ct.find(rule.key) != string::npos) mulContains = rule.multiplier; // last wins
                }
            }
            mul = mulExact.has_value() ? mulExact.value() : (mulContains.has_value() ? mulContains.value() : 1.0);
        }

        const double E = _energy * mul;
        const double F = _force * mul;
        const double Viso = _virialIso * mul;
        const double Vaniso = _virialAniso * mul;

        if (!structure.energySigmaInverse.has_value()) {
            structure.energySigmaInverse = 1.0 / (E * pow(structure.size(), 0.5));
        }

        const double dF = 1.0 / F;
        if (!structure.forceSigmasInverse.has_value()) {
            structure.forceSigmasInverse = vector(structure.size(),  Vector3{dF, dF, dF});
        }

        const double dV_iso = 1.0 / (Viso * pow(structure.size(), 0.5));
        const double dV_aniso = 1.0 / (Vaniso * pow(structure.size(), 0.5));
        if (!structure.virialSigmasInverse.has_value()) {
            // Fill respecting usage in QRGapFit: (0,0), (0,1), (0,2), (1,1), (1,2), (2,2)
            Vector3 row0{dV_iso, dV_aniso, dV_aniso}; // xx, xy, xz
            Vector3 row1{dV_aniso, dV_iso, dV_aniso}; // yx (unused), yy, yz
            Vector3 row2{dV_aniso, dV_aniso, dV_iso}; // zx (unused), zy (unused), zz
            structure.virialSigmasInverse = { row0, row1, row2 };
        }
    }
}
