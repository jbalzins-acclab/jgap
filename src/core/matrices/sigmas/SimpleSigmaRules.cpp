//
// Created by Jegors Balzins on 23.6.2025.
//

#include "core/matrices/sigmas/SimpleSigmaRules.hpp"

namespace jgap {

    SimpleSigmaRules::SimpleSigmaRules(const nlohmann::json &params) {
        _defaultEPerAtom = params["E_per_root_n_atoms"];
        _defaultF = params["F_component"];
        _liquidMultiplier = params["liquid"];
        _shortRangeMultiplier = params["short_range"];
    }

    SimpleSigmaRules::SimpleSigmaRules(double defaultEPerAtom,
                                       double defaultF,
                                       double liquidMultiplier,
                                       double shortRangeMultiplier) {
        _defaultEPerAtom = defaultEPerAtom;
        _defaultF = defaultF;
        _liquidMultiplier = liquidMultiplier;
        _shortRangeMultiplier = shortRangeMultiplier;
    }

    void SimpleSigmaRules::fillSigmas(AtomicStructure &structure) {
        double multiplier = 1.0;
        auto ct = structure.configType.value_or("default");

        if (ct == "isolated_atom") {
            multiplier = 0.001;
        }

        if (ct.contains("liquid") || ct.contains("melt")) {
            multiplier = _liquidMultiplier;
        }

        if (ct.contains("short") || ct.contains("traj") || ct.contains("dimer") || ct.contains("trimer")) {
            multiplier = _shortRangeMultiplier;
        }

        if (!structure.energySigma.has_value()) {
            structure.energySigma = multiplier * _defaultEPerAtom * pow(structure.atoms.size(), 0.5);
            //structure.energySigma = 3.0;
        }

        double dF = multiplier * _defaultF;
        for (auto& atom: structure.atoms) {
            if (!atom.forceSigmas.has_value()) {
                atom.forceSigmas = Vector3{dF, dF, dF};
            } else {
                CurrentLogger::get()->info("skibidi");
            }
        }
    }
}
