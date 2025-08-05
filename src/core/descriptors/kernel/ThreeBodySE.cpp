#include "core/descriptors/kernels/ThreeBodySE.hpp"

#include "utils/Utils.hpp"

namespace jgap {
    ThreeBodySE::ThreeBodySE(double energyScale, double lengthScale) {
        _lengthScale = lengthScale;
        _energyScaleSquared = energyScale * energyScale;
        _inverse2ThetaSq = 1.0 / (2.0 * _lengthScale * _lengthScale);
        _inverseThetaSq = 1.0 / (_lengthScale * _lengthScale);
    }

    ThreeBodySE::ThreeBodySE(const nlohmann::json &params) {
        _lengthScale = params["length_scale"].get<double>();
        _energyScaleSquared = pow(params["energy_scale"].get<double>(), 2);
        _inverse2ThetaSq = 1.0 / (2.0 * _lengthScale * _lengthScale);
        _inverseThetaSq = 1.0 / (_lengthScale * _lengthScale);
    }

    nlohmann::json ThreeBodySE::serialize() {
        return {
            {"length_scale", _lengthScale},
            {"energy_scale", sqrt(_energyScaleSquared)}
        };
    }

    Covariance ThreeBodySE::covariance(const AtomicStructure &structure,
                                       const ThreeBodyKernelIndex &indexes,
                                       const ThreeBodyDescriptorData &sparsePoint) {
        double energy = 0;
        vector forces(structure.size(), Vector3{0, 0, 0});
        array<Vector3, 3> virials{};

        for (ThreeBodyKernelIndexEntity index: indexes) {
            // ---------------------- ENERGY --------------------------------
            energy += covarianceNoCutoffs(index.q, sparsePoint.q)
                        * index.fCut01 * index.fCut02 * sparsePoint.fCut
                        * 2.0/*q_ijk + q_ikj*/;

            // ---------------------- FORCES --------------------------------
            // e.g. dU/dq
            Vector3 gradUWrtQ = gradient(index.q, sparsePoint.q);

            // dU/d{r01, r02, r12} = dU/dq * dq/dr_ij
            Vector3 gradUWrtDistances = Vector3{
                gradUWrtQ.dot(index.dq_k_dr_ij[0]),
                gradUWrtQ.dot(index.dq_k_dr_ij[1]),
                gradUWrtQ.dot(index.dq_k_dr_ij[2])
            } * index.fCut01 * index.fCut02;

            // Product rule with cutoffs
            if (index.fCut01 < 1.0) {
                gradUWrtDistances.x += covarianceNoCutoffs(index.q, sparsePoint.q) * index.dfCut_dr_01 * index.fCut02;
            }
            if (index.fCut02 < 1.0) {
                gradUWrtDistances.y += covarianceNoCutoffs(index.q, sparsePoint.q) * index.fCut01 * index.dfCut_dr_02;
            }
            gradUWrtDistances *= 2.0/*q_ijk + q_ikj */ *  sparsePoint.fCut;

            // chain rule x2 ( remember: gradWrtDistances = d/d{r01, r02, r12} )

            const Vector3 f10 = index.grad_rij_wrt_rj[0] * gradUWrtDistances.x;
            forces[index.atomIndex[0]] += f10;
            forces[index.atomIndex[1]] -= f10;
            virials[0] -= f10 * index.r_ij[0].x;
            virials[1] -= f10 * index.r_ij[0].y;
            virials[2] -= f10 * index.r_ij[0].z;

            const Vector3 f20 = index.grad_rij_wrt_rj[1] * gradUWrtDistances.y;
            forces[index.atomIndex[0]] += f20;
            forces[index.atomIndex[2]] -= f20;
            virials[0] -= f20 * index.r_ij[1].x;
            virials[1] -= f20 * index.r_ij[1].y;
            virials[2] -= f20 * index.r_ij[1].z;

            const Vector3 f21 = index.grad_rij_wrt_rj[2] * gradUWrtDistances.z;
            forces[index.atomIndex[1]] += f21;
            forces[index.atomIndex[2]] -= f21;
            virials[0] -= f21 * index.r_ij[2].x;
            virials[1] -= f21 * index.r_ij[2].y;
            virials[2] -= f21 * index.r_ij[2].z;
        }

        return {energy, forces, virials};
    }

    double ThreeBodySE::covariance(const ThreeBodyDescriptorData &t1,
                                   const ThreeBodyDescriptorData &t2) {
        return covarianceNoCutoffs(t1.q, t2.q) * t1.fCut * t2.fCut;
    }

    double ThreeBodySE::covarianceNoCutoffs(const Vector3 &t1, const Vector3 &t2) const {
        return _energyScaleSquared * exp(-(t1 - t2).square() * _inverse2ThetaSq);
    }

    Vector3 ThreeBodySE::gradient(const Vector3 &changingQ, const Vector3 &constQ) const {
        return (constQ - changingQ) * _inverseThetaSq * covarianceNoCutoffs(changingQ, constQ);
    }
}
