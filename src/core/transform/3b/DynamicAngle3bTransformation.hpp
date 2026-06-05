#ifndef JGAP_DYNAMICANGLE3BTRANSFORMATION_HPP
#define JGAP_DYNAMICANGLE3BTRANSFORMATION_HPP

#include "core/cutoff/CutoffFunction.hpp"
#include "core/transform/ClusterTransformation.hpp"

namespace jgap {

    class DynamicAngle3bTransformation final : public ClusterTransformation<4, 3> {
    public:
        DynamicAngle3bTransformation(std::unique_ptr<CutoffFunction> cutoff)
            : cutoff(std::move(cutoff)) {}

        Descriptor<4> evaluate(const Cluster<3>& triplet) const override {
            Real r01 = triplet.between(0, 1).magnitude;
            Real r02 = triplet.between(0, 2).magnitude;
            Real r12 = triplet.between(1, 2).magnitude;

            Real f_cut_01 = cutoff->evaluate(r01);
            Real f_cut_02 = cutoff->evaluate(r02);

            return {
                .value = {
                    r01 + r02,
                    (r01 - r02) * (r01 - r02),
                    r12,
                    f_cut_01 * f_cut_02
                }
            };
        }

        DescriptorAndDerivatives<4, 3> evaluateAndDifferentiate(const Cluster<3>& triplet) const override {
            Real r01 = triplet.between(0, 1).magnitude;
            Real r02 = triplet.between(0, 2).magnitude;
            Real r12 = triplet.between(1, 2).magnitude;

            auto [f_cut_01, df_cut_01] = cutoff->evaluateAndDifferentiate(r01);
            auto [f_cut_02, df_cut_02] = cutoff->evaluateAndDifferentiate(r02);

            return {
                .value = {
                    r01 + r02,
                    (r01 - r02) * (r01 - r02),
                    r12,
                    f_cut_01 * f_cut_02
                },
                .derivatives = {
                    std::array<Real, 4>{ // wrt r_01
                        1.0,
                        2.0 * (r01 - r02),
                        0.0,
                        df_cut_01 * f_cut_02
                    },
                    std::array<Real, 4>{ // wrt r_02
                        1.0,
                        2.0 * (r02 - r01),
                        0.0,
                        df_cut_02 * f_cut_01
                    },
                    std::array<Real, 4>{ // wrt r_12
                        0.0,
                        0.0,
                        1.0,
                        0.0
                    }
                }
            };
        }

        [[nodiscard]] Cutoffs getCutoffs() const override {
            return Cutoffs{ {3, cutoff->getCutoff()} };
        }

    private:
        std::unique_ptr<CutoffFunction> cutoff;
    };
}

#endif //JGAP_DYNAMICANGLE3BTRANSFORMATION_HPP