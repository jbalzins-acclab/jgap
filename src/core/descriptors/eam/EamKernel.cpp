#include "core/descriptors/eam/EamKernel.hpp"

namespace jgap {

    Covariance EamKernel::covariance(const AtomicStructure &structure,
                                     const EamKernelIndexPerSpecies &indexes) {

        double energy = 0;
        std::vector forces(structure.size(), Vector3{0.0, 0.0, 0.0});
        array<Vector3, 3> virials{};

        for (const auto &index : indexes) {
            energy += valueInternal(index.density);

            const double dU_drho_i = derivativeInternal(index.density);
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

    double EamKernel::value(const double &density) {
        return valueInternal(density);
    }
}