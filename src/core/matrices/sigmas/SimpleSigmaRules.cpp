//
// Created by Jegors Balzins on 23.6.2025.
//

#include "core/matrices/sigmas/SimpleSigmaRules.hpp"

namespace jgap {

    void SimpleSigmaRules::fillSigmas(AtomicStructure &structure) {
        double multiplier = 1.0;
        if (structure.configType.value_or("default").contains("liquid")) {
            multiplier = _liquidMultiplier;
        }
        if (structure.configType.value_or("default").contains("short")) {
            multiplier = _shortRangeMultiplier;
        }

        if (!structure.energySigma.has_value()) {
            structure.energySigma = multiplier * pow(_defaultEPerAtom,1) * pow(structure.atoms.size(), 0.5);
            //structure.energySigma = 3.0;
        }

        double dF = multiplier * _defaultF;
        for (auto& atom: structure.atoms) {
            if (!atom.forceSigmas.has_value()) {
                atom.forceSigmas = Vector3{dF, dF, dF};
            } else {
                Logger::logger->info("skibidi");
            }
        }
    }
}
