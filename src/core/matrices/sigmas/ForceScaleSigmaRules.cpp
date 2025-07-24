#include "core/matrices/sigmas/ForceScaleSigmaRules.hpp"

namespace jgap {
    ForceScaleSigmaRules::ForceScaleSigmaRules(const nlohmann::json &params) {
        // ~ 1.3xSimple:
        _sigmaFmin = params.value("sigma_f_min", 0.0005);
        _referenceF = params.value("reference_f", 0.8);
        _referenceSigmaF = params.value("reference_sigma_f", 0.05);
        _scalingFactor = params.value("scaling_factor", 2);
        _componentWise = params.value("component_wise", false);
        _EFratio = params.value("sigma_E_F_ratio", 1.0 / 50.0);
        _VFratio = params.value("sigma_V_F_ratio", 2.0);
        // TODO: isolated_atom: E:0.008:0.000 F:0.000:0.000
        // Vary ref and scaling
    }

    void ForceScaleSigmaRules::fillSigmas(AtomicStructure &structure) {

        structure.forceSigmasInverse.reset();
        structure.forceSigmasInverse = vector<Vector3>();

        if (structure.configType.value_or("default") == "isolated_atom") {
            structure.energySigmaInverse = 1.0 / 0.00001;
            structure.forceSigmasInverse = {Vector3(1.0 / 0.0005, 1.0 / 0.0005, 1.0 / 0.0005)};
            // TODO: isolated atoms have non-zero virials??
            structure.virials = {};
        }

        double varianceE = 0;
        double sigmaV = 0;
        for (size_t i = 0; i < structure.size(); i++) {
            if (_componentWise) {

                Vector3 sigmaF = {
                    calculateSigma((*structure.forces)[i].x),
                    calculateSigma((*structure.forces)[i].y),
                    calculateSigma((*structure.forces)[i].z)
                };
                structure.forceSigmasInverse->push_back(Vector3{1.0 / sigmaF.x, 1.0 / sigmaF.y, 1.0 / sigmaF.z});

                // maxSigmaF = max(maxSigmaF, sqrt(sigmaX * sigmaX + sigmaY * sigmaY + sigmaZ * sigmaZ));
                varianceE += sigmaF.square() * _EFratio * _EFratio;
                sigmaV += sigmaF.len() * _VFratio;

            } else {
                const double sigmaF = calculateSigma((*structure.forces)[i].len());
                structure.forceSigmasInverse->push_back(Vector3{1.0 / sigmaF, 1.0 / sigmaF, 1.0 / sigmaF});

                // maxSigmaF = max(maxSigmaF, sigmaF);
                varianceE += sigmaF * sigmaF * _EFratio * _EFratio;
                sigmaV += sigmaF * _VFratio;
            }
        }
        sigmaV /= static_cast<double>(structure.size());

        structure.energySigmaInverse = 1.0 / sqrt(varianceE);

        CurrentLogger::get()->debug(format(
            "Sigmas:E{}",
            sqrt(varianceE/static_cast<double>(structure.size()))
            ));

        const double inverseSigmaV = 1.0 / sigmaV;
        structure.virialSigmasInverse = {
            Vector3{inverseSigmaV, inverseSigmaV, inverseSigmaV},
            Vector3{inverseSigmaV, inverseSigmaV, inverseSigmaV},
            Vector3{inverseSigmaV, inverseSigmaV, inverseSigmaV}
        };
    }

    double ForceScaleSigmaRules::calculateSigma(const double F) const {
        return max(_sigmaFmin, F/_referenceF * _scalingFactor * _referenceSigmaF);
    }
}