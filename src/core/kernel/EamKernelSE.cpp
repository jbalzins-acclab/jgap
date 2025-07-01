
#include "core/kernels/EamKernelSE.hpp"

#include <ranges>

namespace jgap {

    double EamKernelSE::covariance(const AtomicStructure &structure,
                                   const EamKernelIndex &indexes,
                                   const double &sparseDensity) {

        double result = 0;

        for (const auto &[density, densityDerivatives] : indexes) {
            result += covariance(density, sparseDensity);
        }

        return result;
    }

    vector<Vector3> EamKernelSE::derivatives(const AtomicStructure &structure,
                                             const EamKernelIndex &indexes,
                                             const double &sparseDensity) {

        vector<Vector3> result(structure.atoms.size());

        for (size_t i = 0; i < indexes.size(); i++) {
            double dK_drho_i = derivative(indexes[i].density, sparseDensity);
            auto atom = structure.atoms[i];

            for (auto &[neighbourData, d_rho_i_dr_ij]: indexes[i].densityDerivatives) {
                Vector3 displacement = structure.atoms[neighbourData.index].position + neighbourData.offset
                                    - atom.position;
                Vector3 df = displacement.normalize() * d_rho_i_dr_ij * dK_drho_i;
                result[i] = result[i] - df;
                result[neighbourData.index] = result[neighbourData.index] + df;
            }
        }

        return result;
    }

    double EamKernelSE::covariance(const double &density1, const double &density2) {
        return _energyScaleSquared * exp(-pow(density1-density2, 2) * _inverse2ThetaSq);
    }

    double EamKernelSE::derivative(const double &changingDensity, const double &constantDensity) {
        return (constantDensity - changingDensity) * 2/*compensate constant*/
                    * _inverse2ThetaSq * covariance(changingDensity, constantDensity);
    }
}
