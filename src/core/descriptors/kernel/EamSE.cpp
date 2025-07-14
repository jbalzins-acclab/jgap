
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

        for (const auto &[density, densityDerivatives] : indexes) {
            result += covariance(density, sparseDensity);
        }

        return result;
    }

    vector<Vector3> EamSE::derivatives(const AtomicStructure &structure,
                                       const EamKernelIndexPerSpecies &indexes,
                                       const double &sparseDensity) {

        vector<Vector3> result(structure.atoms.size(), {0.0, 0.0, 0.0});

        for (size_t i = 0; i < indexes.size(); i++) {
            double dK_drho_i = derivative(indexes[i].density, sparseDensity);
            CurrentLogger::get()->debug(format("dK_drho_i {}: {} !! {}", dK_drho_i, indexes[i].density, sparseDensity));
            auto atom = structure.atoms[i];

            for (auto &[neighbourData, d_rho_i_dr_ij]: indexes[i].densityDerivatives) {
                Vector3 displacement = structure.atoms[neighbourData.index].position + neighbourData.offset
                                    - atom.position;
                Vector3 df = displacement.normalize() * d_rho_i_dr_ij * dK_drho_i;
            CurrentLogger::get()->debug(format("df {}: {} !! {}", df.toString(), indexes[i].density, sparseDensity));
            CurrentLogger::get()->debug(format("disp {}: {} !! {}", displacement.normalize().toString(), indexes[i].density, sparseDensity));
            CurrentLogger::get()->debug(format("d_rho_i_dr_ij {}: {} !! {}", d_rho_i_dr_ij, indexes[i].density, sparseDensity));
                result[i] = result[i] - df;
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
