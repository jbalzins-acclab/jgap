#include "core/descriptors/3b/ThreeBodyKernel.hpp"

namespace jgap {
    double ThreeBodyKernel::value(const ThreeBodyDescriptorData &t) {
        return valueInternal(t.q) * t.fCut;
    }

    Covariance ThreeBodyKernel::covariance(const AtomicStructure &structure,
                                           const ThreeBodyIndex &index) {
        double energy = 0;
        std::vector forces(structure.size(), Vector3{0, 0, 0});
        array<Vector3, 3> virials{};

        for (auto [atomIndexes,
                   r_ij, grad_rij_wrt_rj,
                   fCut01, fCut02, dfCut_dr_01, dfCut_dr_02,
                   q, dq_k_dr_ij]: index) {
            // ---------------------- ENERGY --------------------------------
            energy += valueInternal(q)
                        * fCut01 * fCut02
                        * 2.0/*q_ijk + q_ikj*/;

            // ---------------------- FORCES --------------------------------
            // e.g. dU/dq
            Vector3 gradUWrtQ = gradientInternal(q);

            // dU/d{r01, r02, r12} = dU/dq * dq/dr_ij
            Vector3 gradUWrtDistances = Vector3{
                gradUWrtQ.dot(dq_k_dr_ij[0]),
                gradUWrtQ.dot(dq_k_dr_ij[1]),
                gradUWrtQ.dot(dq_k_dr_ij[2])
            } * fCut01 * fCut02;

            // Product rule with cutoffs
            if (fCut01 < 1.0) {
                gradUWrtDistances.x += valueInternal(q) * dfCut_dr_01 * fCut02;
            }
            if (fCut02 < 1.0) {
                gradUWrtDistances.y += valueInternal(q) * fCut01 * dfCut_dr_02;
            }
            gradUWrtDistances *= 2.0/*q_ijk + q_ikj */;

            // chain rule x2 ( remember: gradWrtDistances = d/d{r01, r02, r12} )

            const Vector3 f10 = grad_rij_wrt_rj[0] * gradUWrtDistances.x;
            forces[atomIndexes[0]] += f10;
            forces[atomIndexes[1]] -= f10;
            virials[0] -= f10 * r_ij[0].x;
            virials[1] -= f10 * r_ij[0].y;
            virials[2] -= f10 * r_ij[0].z;

            const Vector3 f20 = grad_rij_wrt_rj[1] * gradUWrtDistances.y;
            forces[atomIndexes[0]] += f20;
            forces[atomIndexes[2]] -= f20;
            virials[0] -= f20 * r_ij[1].x;
            virials[1] -= f20 * r_ij[1].y;
            virials[2] -= f20 * r_ij[1].z;

            const Vector3 f21 = grad_rij_wrt_rj[2] * gradUWrtDistances.z;
            forces[atomIndexes[1]] += f21;
            forces[atomIndexes[2]] -= f21;
            virials[0] -= f21 * r_ij[2].x;
            virials[1] -= f21 * r_ij[2].y;
            virials[2] -= f21 * r_ij[2].z;
        }

        return {energy, forces, virials};
    }
}
