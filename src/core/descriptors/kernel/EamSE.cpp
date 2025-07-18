
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

    Covariance EamSE::covariance(const AtomicStructure &structure,
                                 const EamKernelIndexPerSpecies &indexes,
                                 const double &sparseDensity) {

        double energy = 0;
        vector forces(structure.size(), Vector3{0.0, 0.0, 0.0});

        for (const auto &index : indexes) {
            energy += covarianceNoCutoffs(index.density, sparseDensity);

            const double dK_drho_i = derivative(index.density, sparseDensity);
            auto atomPosition = structure.positions[index.atAtomIndex];

            for (auto &[neighbourData, d_rho_i_dr_ij]: index.densityDerivatives) {
                const Vector3 displacement = structure.positions[neighbourData.index] + neighbourData.offset
                                             - atomPosition;
                const Vector3 df = displacement.normalize() * d_rho_i_dr_ij * dK_drho_i;
                forces[index.atAtomIndex] -= df;
                forces[neighbourData.index] += df;
            }
        }

        return {energy, forces};
    }

    double EamSE::covarianceNoCutoffs(const double &density1, const double &density2) const {
        return _energyScaleSquared * exp(-pow(density1-density2, 2) * _inverse2ThetaSq);
    }
    double EamSE::covariance(const double &density1, const double &density2) {
        return covarianceNoCutoffs(density1, density2);
    }

    double EamSE::derivative(const double &changingDensity, const double &constantDensity) const {
        return (constantDensity - changingDensity) * (2/*compensate constant*/* _inverse2ThetaSq)
                     * covarianceNoCutoffs(changingDensity, constantDensity);
    }
}
