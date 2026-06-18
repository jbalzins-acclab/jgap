#ifndef JGAP_ANGLE3BTRANSFORMATION_HPP
#define JGAP_ANGLE3BTRANSFORMATION_HPP

#include "core/cutoff/CutoffFunction.hpp"
#include "core/transform/ClusterTransformation.hpp"

namespace jgap {

    class Angle3bTransformation final : public ClusterTransformation<4, 3> {
    public:

        Angle3bTransformation(const ValuePtr<CutoffFunction>& cutoff)
            : cutoff(cutoff) {}

        Descriptor<4> evaluate(const Cluster<3>& triplet) const override {
            Real r01 = triplet.between(0, 1);
            Real r02 = triplet.between(0, 2);
            Real r12 = triplet.between(1, 2);

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

        NBodyDescriptor<4, 3> evaluateAndDifferentiate(const Cluster<3>& triplet) const override {
            Real r01 = triplet.between(0, 1);
            Real r02 = triplet.between(0, 2);
            Real r12 = triplet.between(1, 2);

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
                    std::array{ // wrt r_01
                        1.0,
                        2.0 * (r01 - r02),
                        0.0,
                        df_cut_01 * f_cut_02
                    },
                    std::array{ // wrt r_02
                        1.0,
                        2.0 * (r02 - r01),
                        0.0,
                        df_cut_02 * f_cut_01
                    },
                    std::array{ // wrt r_12
                        0.0,
                        0.0,
                        1.0,
                        0.0
                    }
                }
            };
        }

        Cutoffs getCutoffs() const override {
            return Cutoffs{ {3, cutoff->getCutoff()} };
        }

        ValuePtr<CutoffFunction> getCutoff() const {
            return cutoff;
        }

        Real symmetryFactor() const override {
            return 2.0;
        }

        std::unique_ptr<ClusterTransformation<4, 3>> clone() const override {
            return std::make_unique<Angle3bTransformation>(*this);
        }

    private:
        ValuePtr<CutoffFunction> cutoff;
    };
}

#endif