//
// Created by Jegors Balzins on 21.7.2025.
//

#include "core/matrices/sigmas/ForceScaleSigmaRules.hpp"

namespace jgap {
    ForceScaleSigmaRules::ForceScaleSigmaRules(const nlohmann::json &params) {
        _Fmin = params.value("f_min", 0.02);
        _scalingFactor = params.value("scaling_factor", 0.1);
        _componentWise = params.value("component_wise", false);
    }

    void ForceScaleSigmaRules::fillSigmas(AtomicStructure &structure) {

        structure.forceSigmasInverse.reset();
        structure.forceSigmasInverse = vector<Vector3>();

        //double maxSigmaF = 0.0;
        double meanSigmaF = 0.0;
        for (size_t i = 0; i < structure.size(); i++) {
            if (_componentWise) {

                const double sigmaX = calculateSigma((*structure.forces)[i].x);
                const double sigmaY = calculateSigma((*structure.forces)[i].y);
                const double sigmaZ = calculateSigma((*structure.forces)[i].z);
                structure.forceSigmasInverse->push_back(Vector3{1.0 / sigmaX, 1.0 / sigmaY, 1.0 / sigmaZ});

                // maxSigmaF = max(maxSigmaF, sqrt(sigmaX * sigmaX + sigmaY * sigmaY + sigmaZ * sigmaZ));
                meanSigmaF += sqrt(sigmaX * sigmaX + sigmaY * sigmaY + sigmaZ * sigmaZ);

            } else {
                const double sigmaF = calculateSigma((*structure.forces)[i].norm());
                structure.forceSigmasInverse->push_back(Vector3{1.0 / sigmaF, 1.0 / sigmaF, 1.0 / sigmaF});

                // maxSigmaF = max(maxSigmaF, sigmaF);
                meanSigmaF += sigmaF;
            }
        }
        meanSigmaF /= static_cast<double>(structure.size());

        const double sigmaE = max(0.002 * sqrt(structure.size()), meanSigmaF * _E_F_RATIO * sqrt(structure.size()));
        structure.energySigmaInverse = 1.0 / sigmaE;

        CurrentLogger::get()->debug(format(
            "Sigmas:F{}-E{}(n={})",
            meanSigmaF, sigmaE, structure.size()
            ));

        const double dV = 1.0 / (meanSigmaF * _V_F_RATIO);
        structure.virialSigmasInverse = {
            Vector3{dV, dV, dV},
            Vector3{dV, dV, dV},
            Vector3{dV, dV, dV}
        };
    }

    double ForceScaleSigmaRules::calculateSigma(const double F) const {
        if (abs(F) < _Fmin) {
            return _scalingFactor * _Fmin;
        }
        return _scalingFactor * abs(F);
    }
}