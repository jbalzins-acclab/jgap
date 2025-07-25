#include "core/matrices/sigmas/SimpleSigmaRules.hpp"

namespace jgap {

    SimpleSigmaRules::SimpleSigmaRules(const nlohmann::json &params) {
        _defaultEPerAtom = params["E_per_root_n_atoms"];
        _defaultF = params.value("F_component", _defaultEPerAtom * 50.0);
        _defaultVirials = params.value("virials", _defaultEPerAtom * 100.0);
        _liquidMultiplier = params["liquid"];
        _shortRangeMultiplier = params["short_range"];
    }

    void SimpleSigmaRules::fillSigmas(AtomicStructure &structure) {
        double multiplier = 1.0;
        const auto ct = structure.configType.value_or("default");

        if (ct == "isolated_atom") {
            multiplier = 0.001;
        }

        if (ct.contains("liquid") || ct.contains("melt")) {
            multiplier = _liquidMultiplier;
        }

        if (ct.contains("short") || ct.contains("traj") || ct.contains("low_volume")
            || ct.contains("dimer") || ct.contains("trimer")) {
            multiplier = _shortRangeMultiplier;
        }

        if (!structure.energySigmaInverse.has_value()) {
            structure.energySigmaInverse = 1.0 / (multiplier * _defaultEPerAtom * pow(structure.size(), 0.5));
        }

        const double dF = 1.0 / (multiplier * _defaultF);
        if (!structure.forceSigmasInverse.has_value()) {
            structure.forceSigmasInverse = vector(structure.size(),  Vector3{dF, dF, dF});
        }

        const double dV = 1.0 / (multiplier * _defaultVirials * pow(structure.size(), 0.5));
        if (!structure.virialSigmasInverse.has_value()) {
            structure.virialSigmasInverse = {
                Vector3{dV, dV, dV},
                Vector3{dV, dV, dV},
                Vector3{dV, dV, dV}
            };
        }
    }
}
