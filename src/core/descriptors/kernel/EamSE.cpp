
#include "core/descriptors/kernels/EamSE.hpp"

#include <ranges>

namespace jgap {
    EamSE::EamSE(double energyScale, double lengthScale) {
        _lengthScale = lengthScale;
        _energyScaleSquared = energyScale * energyScale;
        _inverse2ThetaSq = 1.0 / (2.0 * _lengthScale * _lengthScale);
        _inverseThetaSq = 1.0 / (_lengthScale * _lengthScale);
    }

    EamSE::EamSE(const nlohmann::json &params) {
        _lengthScale = params["length_scale"].get<double>();
        _energyScaleSquared = pow(params["energy_scale"].get<double>(), 2);
        _inverse2ThetaSq = 1.0 / (2.0 * _lengthScale * _lengthScale);
        _inverseThetaSq = 1.0 / (_lengthScale * _lengthScale);
    }

    nlohmann::json EamSE::serialize() {
        return {
            {"length_scale", _lengthScale},
            {"energy_scale", sqrt(_energyScaleSquared)}
        };
    }

    Covariance EamSE::covariance(const AtomicStructure &structure,
                                 const EamKernelIndexPerSpecies &indexes,
                                 const double &sparseDensity) {

        double energy = 0;
        vector forces(structure.size(), Vector3{0.0, 0.0, 0.0});
        array<Vector3, 3> virials{};

        for (const auto &index : indexes) {
            energy += covariance(index.density, sparseDensity);

            const double dU_drho_i = derivative(index.density, sparseDensity);
            auto atomPosition = structure.positions[index.atAtomIndex];

            for (auto &[neighbourData, d_rho_i_dr_ij]: index.densityDerivatives) {
                const Vector3 r01 = structure.positions[neighbourData.index] + neighbourData.offset - atomPosition;
                const Vector3 f10 = r01.normalize() * d_rho_i_dr_ij * dU_drho_i;
                forces[index.atAtomIndex] += f10;
                forces[neighbourData.index] -= f10;

                // x2 since r10.x * f10.x = r01.x * f01.x
                virials[0] += f10 * r01.x;
                virials[1] += f10 * r01.y;
                virials[2] += f10 * r01.z;
            }
        }

        return {energy, forces, virials};
    }

    double EamSE::covariance(const double &density1, const double &density2) {
        return _energyScaleSquared * exp(-pow(density1-density2, 2) * _inverse2ThetaSq);
    }

    double EamSE::derivative(const double &changingDensity, const double &constantDensity) {
        return (constantDensity - changingDensity) * _inverseThetaSq * covariance(changingDensity, constantDensity);
    }
}
