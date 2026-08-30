#ifndef JGAP_MEAMTRANSFORMATION_HPP
#define JGAP_MEAMTRANSFORMATION_HPP

#include "jgap/core/cutoff/CutoffFunction.hpp"
#include "jgap/core/transform/nbody/3b/ThreeBodyTransformation.hpp"

namespace jgap {

    class MeamTransformation final : public ThreeBodyTransformation<3> {
    public:
        MeamTransformation(const ValuePtr<CutoffFunction>& cutoff) : cutoff(cutoff) {}

        Descriptor<3> evaluate(const Cluster3& cluster) const override {
            return ThreeBodyTransformation<3>::evaluate(cluster);
        }

        ThreeBodyDescriptor<3> evaluateAndDifferentiate(const Cluster3& cluster) const override final {
            Real r01 = cluster.separation01.magnitude;
            Real r02 = cluster.separation02.magnitude;
            Real r12 = cluster.separation12.magnitude;

            auto [f1, df1] = cutoff->evaluateAndDifferentiate(r01);
            auto [f2, df2] = cutoff->evaluateAndDifferentiate(r02);

            Real r01_inv = 1.0_r / r01;
            Real r02_inv = 1.0_r / r02;
            Real r01_sq = r01 * r01;
            Real r02_sq = r02 * r02;
            Real r12_sq = r12 * r12;

            Real c = (r01_sq + r02_sq - r12_sq) * 0.5_r * r01_inv * r02_inv;
            // Prevent numerical issues if c goes slightly out of bounds
            if (c > 1.0_r) {
                c = 1.0_r;
            } else if (c < -1.0_r) {
                c = -1.0_r;
            }

            Real dc_dr01 = (r01_sq - r02_sq + r12_sq) * 0.5_r * r01_inv * r01_inv * r02_inv;
            Real dc_dr02 = (r02_sq - r01_sq + r12_sq) * 0.5_r * r02_inv * r02_inv * r01_inv;
            Real dc_dr12 = -r12 * r01_inv * r02_inv;

            Real c2 = c * c;
            Real c3 = c2 * c;

            Real L1 = c;
            Real L2 = c2 - 1.0_r / 3.0_r;
            Real L3 = c3 - 0.6_r * c;

            Real dL1_dc = 1.0_r;
            Real dL2_dc = 2.0_r * c;
            Real dL3_dc = 3.0_r * c2 - 0.6_r;

            Real f1f2 = f1 * f2;

            const auto& dir01 = cluster.separation01.direction;
            const auto& dir02 = cluster.separation02.direction;
            const auto& dir12 = cluster.separation12.direction;

            ThreeBodyDescriptor<3> desc;
            desc.value[0] = f1f2 * L1;
            desc.value[1] = f1f2 * L2;
            desc.value[2] = f1f2 * L3;

            Real df1f2 = df1 * f2;
            Real f1df2 = f1 * df2;

            std::array<Real, 3> dq_dr01 = {
                df1f2 * L1 + f1f2 * dL1_dc * dc_dr01,
                df1f2 * L2 + f1f2 * dL2_dc * dc_dr01,
                df1f2 * L3 + f1f2 * dL3_dc * dc_dr01
            };

            std::array<Real, 3> dq_dr02 = {
                f1df2 * L1 + f1f2 * dL1_dc * dc_dr02,
                f1df2 * L2 + f1f2 * dL2_dc * dc_dr02,
                f1df2 * L3 + f1f2 * dL3_dc * dc_dr02
            };

            std::array<Real, 3> dq_dr12 = {
                f1f2 * dL1_dc * dc_dr12,
                f1f2 * dL2_dc * dc_dr12,
                f1f2 * dL3_dc * dc_dr12,
            };

            for (size_t d = 0; d < 3; ++d) {
                desc.grad_r1[d] = dq_dr01[d] * dir01 - dq_dr12[d] * dir12;
                desc.grad_r2[d] = dq_dr02[d] * dir02 + dq_dr12[d] * dir12;
            }

            return desc;
        }

        Cutoffs getCutoffs() const override { return Cutoffs{{3, cutoff->getCutoff()}}; }
        bool isRotationallyInvariant() const override { return true; }
        bool isSwapInvariant(size_t idx1, size_t idx2) const override {
            return (idx1 == 1 && idx2 == 2) || (idx1 == 2 && idx2 == 1);
        }

        ValuePtr<CutoffFunction> getCutoffFunction() const { return cutoff; }

        MeamTransformation* clone() const override { return new MeamTransformation(*this); }

    private:
        ValuePtr<CutoffFunction> cutoff;
    };
}

#endif
