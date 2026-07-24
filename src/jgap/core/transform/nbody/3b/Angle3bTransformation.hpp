#ifndef JGAP_ANGLE3BTRANSFORMATION_HPP
#define JGAP_ANGLE3BTRANSFORMATION_HPP

#include "ThreeBodyTransformation.hpp"
#include "jgap/core/cutoff/CutoffFunction.hpp"

namespace jgap {

    class Angle3bTransformation final : public ThreeBodyTransformation<4> {
    public:
        Angle3bTransformation(const ValuePtr<CutoffFunction>& cutoff) : cutoff(cutoff) {}

        Descriptor<4> evaluate(const Cluster3& triplet) const override {
            Real r01 = triplet.separationBetween(0, 1);
            Real r02 = triplet.separationBetween(0, 2);
            Real r12 = triplet.separationBetween(1, 2);

            Real f_cut_01 = cutoff->evaluate(r01);
            Real f_cut_02 = cutoff->evaluate(r02);

            return {r01 + r02, (r01 - r02) * (r01 - r02), r12, f_cut_01 * f_cut_02};
        }

        ThreeBodyDescriptor<4> evaluateAndDifferentiate(const Cluster3& triplet) const override {
            Real r01 = triplet.separationBetween(0, 1);
            Real r02 = triplet.separationBetween(0, 2);
            Real r12 = triplet.separationBetween(1, 2);

            auto [f_cut_01, df_cut_01] = cutoff->evaluateAndDifferentiate(r01);
            auto [f_cut_02, df_cut_02] = cutoff->evaluateAndDifferentiate(r02);

            // clang-format off
            return {
                .value = {
                    r01 + r02,
                    (r01 - r02) * (r01 - r02),
                    r12,
                    f_cut_01 * f_cut_02,
                },
                .derivatives = {
                    std::array{
                        // wrt r_01
                        1.0,
                        2.0 * (r01 - r02),
                        0.0,
                        df_cut_01 * f_cut_02,
                    },
                    std::array{
                        // wrt r_02
                        1.0,
                        2.0 * (r02 - r01),
                        0.0,
                        df_cut_02 * f_cut_01,
                    },
                    std::array{
                        // wrt r_12
                        0.0,
                        0.0,
                        1.0,
                        0.0,
                    },
                }
            };
            // clang-format on
        }

        Cutoffs getCutoffs() const override { return Cutoffs{{3, cutoff->getCutoff()}}; }

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
