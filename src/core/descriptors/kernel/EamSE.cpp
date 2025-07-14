
#include "core/descriptors/kernels/EamSE.hpp"

#include <ranges>

namespace jgap {
    EamSE::EamSE(const nlohmann::json &params) {
        _lengthScale = params["length_scale"].get<double>();
        _energyScaleSquared = pow(params["energy_scale"].get<double>(), 2);
        _inverse2ThetaSq = 1.0 / (2.0 * _lengthScale * _lengthScale);
    }

    nlohmann::json EamSE::serialize() {
        return {
            {"length_scale", _lengthScale},
            {"energy_scale", sqrt(_energyScaleSquared)}
        };
    }

    EamSE::EamSE(const double energyScale, const double lengthScale): _lengthScale(lengthScale) {
        _inverse2ThetaSq = 1.0 / (2.0 * _lengthScale * _lengthScale);
        _energyScaleSquared = energyScale * energyScale;
    }

    double EamSE::covariance(const AtomicStructure &structure,
                             const EamKernelIndexPerSpecies &indexes,
                             const double &sparseDensity) {

        double result = 0;

        for (const auto &[atAtomIndex, density, densityDerivatives] : indexes) {
            result += covariance(density, sparseDensity);
        }

        return result;
    }

    vector<Vector3> EamSE::derivatives(const AtomicStructure &structure,
                                       const EamKernelIndexPerSpecies &indexes,
                                       const double &sparseDensity) {

        vector<Vector3> result(structure.atoms.size(), {0.0, 0.0, 0.0});

        for (const auto &index : indexes) {
            const double dK_drho_i = derivative(index.density, sparseDensity);
            auto atom = structure.atoms[index.atAtomIndex];

            for (auto &[neighbourData, d_rho_i_dr_ij]: index.densityDerivatives) {
                Vector3 displacement = structure.atoms[neighbourData.index].position + neighbourData.offset
                                        - atom.position;
                Vector3 df = displacement.normalize() * d_rho_i_dr_ij * dK_drho_i;
                result[index.atAtomIndex] = result[index.atAtomIndex] - df;
                result[neighbourData.index] = result[neighbourData.index] + df;
            }
        }

        return result;
    }

    double EamSE::covariance(const double &density1, const double &density2) {
        return _energyScaleSquared * exp(-pow(density1-density2, 2) * _inverse2ThetaSq);
    }

    double EamSE::derivative(const double &changingDensity, const double &constantDensity) {
        return (constantDensity - changingDensity) * 2/*compensate constant*/
                    * _inverse2ThetaSq * covariance(changingDensity, constantDensity);
    }
}
