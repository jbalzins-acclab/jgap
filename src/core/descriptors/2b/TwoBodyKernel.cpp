#include "core/descriptors/2b/TwoBodyKernel.hpp"

namespace jgap {
    Covariance TwoBodyKernel::covariance(const AtomicStructure &structure, const TwoBodyIndex &index) {
        double energy = 0;
        vector<Vector3> forces(structure.size(), {0, 0, 0});
        array<Vector3, 3> virials{};

        for (const auto &[atomIndex0, atomIndex1, r01, r, fCut, dCut_dr]: index) {
            // ---------------------- ENERGY --------------------------------
            auto dE = valueInternal(r) * fCut;
            if (atomIndex0 != atomIndex1) dE *= 2.0; // K(r_ij,)+K(r_ji)(?)
            // else *= 1.0 --- since both r_ii(+offset) and r_ii(-offset) are in index
            energy += dE;

            // ---------------------- FORCES --------------------------------
            double dE_dr = derivativeInternal(r) * fCut;

            if (fCut < 1.0) {
                dE_dr += valueInternal(r) * dCut_dr;
            }

            auto f10 = r01.normalize() * dE_dr;
            if (atomIndex0 != atomIndex1) {
                f10 *= 2.0;
                forces[atomIndex0] += f10;
                forces[atomIndex1] -= f10;
            }

            virials[0] -= f10 * r01.x;
            virials[1] -= f10 * r01.y;
            virials[2] -= f10 * r01.z;
        }

        return {energy, forces, virials};
    }

    double TwoBodyKernel::value(const TwoBodyDescriptorData &r) {
        return valueInternal(r.r) * r.fCut;
    }
}
