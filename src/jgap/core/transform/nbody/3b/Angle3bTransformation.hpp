#ifndef JGAP_ANGLE3BTRANSFORMATION_HPP
#define JGAP_ANGLE3BTRANSFORMATION_HPP

#include "ThreeBodyTransformation.hpp"
#include "jgap/core/cutoff/CutoffFunction.hpp"

namespace jgap {

    class Angle3bTransformation final : public ThreeBodyTransformation<4> {
    public:
        Angle3bTransformation(const ValuePtr<CutoffFunction>& cutoff) : cutoff(cutoff) {}

        Descriptor<4> evaluate(const Cluster3& triplet) const override {
            return ThreeBodyTransformation<4>::evaluate(triplet);
        }

        ThreeBodyDescriptor<4> evaluateAndDifferentiate(const Cluster3& triplet) const override final {
            Real r01 = triplet.separation01.magnitude;
            Real r02 = triplet.separation02.magnitude;
            Real r12 = triplet.separation12.magnitude;

            const auto& dir01 = triplet.separation01.direction;
            const auto& dir02 = triplet.separation02.direction;
            const auto& dir12 = triplet.separation12.direction;

            auto [f_cut_01, df_cut_01] = cutoff->evaluateAndDifferentiate(r01);
            auto [f_cut_02, df_cut_02] = cutoff->evaluateAndDifferentiate(r02);

            return {
                .value = {
                    r01 + r02,
                    (r01 - r02) * (r01 - r02),
                    r12,
                    f_cut_01 * f_cut_02,
                },
                .grad_r1 = {
                    dir01,
                    2.0_r * (r01 - r02) * dir01,
                    -dir12,
                    (df_cut_01 * f_cut_02) * dir01,
                },
                .grad_r2 = {
                    dir02,
                    2.0_r * (r02 - r01) * dir02,
                    dir12,
                    (df_cut_02 * f_cut_01) * dir02,
                }
            };
        }

        Cutoffs getCutoffs() const override { return Cutoffs{ { 3, cutoff->getCutoff() } }; }
        bool isRotationallyInvariant() const override { return true; }

        ValuePtr<CutoffFunction> getCutoffFunction() const { return cutoff; }

        Angle3bTransformation* clone() const override { return new Angle3bTransformation(*this); }

        bool isSwapInvariant(size_t idx1, size_t idx2) const override {
            return (idx1 == 1 && idx2 == 2) || (idx1 == 2 && idx2 == 1);
        }

    private:
        ValuePtr<CutoffFunction> cutoff;
    };
}

#endif
