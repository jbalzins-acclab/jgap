#ifndef JGAP_CUTOFFJK3BTRANSFORMATION_HPP
#define JGAP_CUTOFFJK3BTRANSFORMATION_HPP

#include <algorithm>
#include <optional>

#include "jgap/core/cutoff/CutoffFunction.hpp"
#include "jgap/core/transform/nbody/3b/ThreeBodyTransformation.hpp"

namespace jgap {

    class CutoffJK3bTransformation final : public ThreeBodyTransformation<4> {
    public:
        CutoffJK3bTransformation(
            const ValuePtr<CutoffFunction>& main_cutoff,
            const std::optional<ValuePtr<CutoffFunction>>& cutoff_12_opt = std::nullopt
        ) :
            main_cutoff(main_cutoff), cutoff_12(cutoff_12_opt.value_or(main_cutoff)) {}

        Descriptor<4> evaluate(const Cluster3& triplet) const override {
            return ThreeBodyTransformation<4>::evaluate(triplet);
        }

        ThreeBodyDescriptor<4> evaluateAndDifferentiate(const Cluster3& triplet) const override final {
            Real r01 = triplet.separation01().magnitude;
            Real r02 = triplet.separation02().magnitude;
            Real r12 = triplet.separation12().magnitude;

            auto [f_cut_01, df_cut_01] = main_cutoff->evaluateAndDifferentiate(r01);
            auto [f_cut_02, df_cut_02] = main_cutoff->evaluateAndDifferentiate(r02);
            auto [f_cut_12, df_cut_12] = cutoff_12->evaluateAndDifferentiate(r12);

            return {
                .value =
                    {
                        r01 + r02,
                        (r01 - r02) * (r01 - r02),
                        r12,
                        f_cut_01 * f_cut_02 * f_cut_12,
                    },
                .derivatives = {
                    std::array{
                        // wrt r_01
                        1.0_r,
                        2.0_r * (r01 - r02),
                        0.0_r,
                        df_cut_01 * f_cut_02 * f_cut_12,
                    },
                    std::array{
                        // wrt r_02
                        1.0_r,
                        2.0_r * (r02 - r01),
                        0.0_r,
                        df_cut_02 * f_cut_01 * f_cut_12,
                    },
                    std::array{
                        // wrt r_12
                        0.0_r,
                        0.0_r,
                        1.0_r,
                        df_cut_12 * f_cut_01 * f_cut_02,
                    },
                }
            };
        }

        Cutoffs getCutoffs() const override {
            return Cutoffs{{3, std::max(main_cutoff->getCutoff(), cutoff_12->getCutoff())}};
        }

        const ValuePtr<CutoffFunction>& getMainCutoffFunction() const { return main_cutoff; }
        const ValuePtr<CutoffFunction>& getCutoff12Function() const { return cutoff_12; }

        CutoffJK3bTransformation* clone() const override { return new CutoffJK3bTransformation(*this); }

        bool isSwapInvariant(size_t idx1, size_t idx2) const override {
            return (idx1 == 1 && idx2 == 2) || (idx1 == 2 && idx2 == 1);
        }

    private:
        ValuePtr<CutoffFunction> main_cutoff;
        ValuePtr<CutoffFunction> cutoff_12;
    };
}

#endif
