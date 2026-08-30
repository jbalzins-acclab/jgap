#ifndef JGAP_DISTANCES3B_HPP
#define JGAP_DISTANCES3B_HPP

#include "../../../../core/transform/nbody/3b/ThreeBodyTransformation.hpp"
#include "jgap/core/cutoff/CutoffFunction.hpp"

namespace jgap {

    class Distances3bTransformation final : public ThreeBodyTransformation<4> {
    public:
        Distances3bTransformation(const ValuePtr<CutoffFunction>& cutoff) : cutoff(cutoff) {}

        Descriptor<4> evaluate(const Cluster3& triplet) const override {
            return ThreeBodyTransformation<4>::evaluate(triplet);
        }

        ThreeBodyDescriptor<4> evaluateAndDifferentiate(const Cluster3& triplet) const override final {
            Real r01 = triplet.separation01.magnitude;
            Real r02 = triplet.separation02.magnitude;
            Real r12 = triplet.separation12.magnitude;

            auto [f_cut_01, df_cut_01] = cutoff->evaluateAndDifferentiate(r01);
            auto [f_cut_02, df_cut_02] = cutoff->evaluateAndDifferentiate(r02);

            const auto& dir01 = triplet.separation01.direction;
            const auto& dir02 = triplet.separation02.direction;
            const auto& dir12 = triplet.separation12.direction;

            return {
                .value = {
                    r01,
                    r02,
                    r12,
                    f_cut_01 * f_cut_02,
                },
                .grad_r1 = {
                    dir01,
                    Vector3{},
                    -dir12,
                    (df_cut_01 * f_cut_02) * dir01,
                },
                .grad_r2 = {
                    Vector3{},
                    dir02,
                    dir12,
                    (df_cut_02 * f_cut_01) * dir02,
                }
            };
        }

        Cutoffs getCutoffs() const override { return Cutoffs{ { 3, cutoff->getCutoff() } }; }
        bool isRotationallyInvariant() const override { return true; }

        ValuePtr<CutoffFunction> getCutoffFunction() const { return cutoff; }

        Distances3bTransformation* clone() const override { return new Distances3bTransformation(*this); }

        bool isSwapInvariant(size_t idx1, size_t idx2) const override { return false; }

    private:
        ValuePtr<CutoffFunction> cutoff;
    };
}

#endif
