#ifndef JGAP_ANGLE3BTRANSFORMER_HPP
#define JGAP_ANGLE3BTRANSFORMER_HPP

#include "core/cutoff/CutoffFunction.hpp"
#include "../ClusterTransformation.hpp"
#include "core/cutoff/CosCutoff.hpp"

namespace jgap {

    template<typename TCutoff = CosCutoff>
    requires std::derived_from<TCutoff, CutoffFunction>
    class Angle3bTransformation : public ClusterTransformation<4, 3> {
    public:
        static constexpr size_t Dim = 4;
        static constexpr size_t ClusterSize = 3;

        Angle3bTransformation() = default;
        Angle3bTransformation(TCutoff cutoff) : cutoff(cutoff) { }

        Descriptor<4> evaluate(const Cluster<3> &triplet) const override;

        DescriptorAndDerivatives<4, 3> evaluateAndDifferentiate(const Cluster<3> &triplet) const override;

        Cutoffs getCutoffs() const override{
            return Cutoffs{ {3, cutoff.getCutoff() } };
        }

    private:
        TCutoff cutoff;
    };

    static_assert(CClusterTransformation<Angle3bTransformation<>>);

    template<typename TCutoff> requires std::derived_from<TCutoff, CutoffFunction>
    Descriptor<4> Angle3bTransformation<TCutoff>::evaluate(const Cluster<3> &triplet) const {
        Real r01 = triplet.between(0, 1).magnitude;
        Real r02 = triplet.between(0, 2).magnitude;
        Real r12 = triplet.between(1, 2).magnitude;

        Real f_cut_01 = cutoff.evaluate(r01);
        Real f_cut_02 = cutoff.evaluate(r02);

        return {
            .value = {
                r01 + r02,
                (r01 - r02) * (r01 - r02),
                r12,
                f_cut_01 * f_cut_02
            }
        };
    }

    template<typename TCutoff> requires std::derived_from<TCutoff, CutoffFunction>
    DescriptorAndDerivatives<4, 3> Angle3bTransformation<TCutoff>::evaluateAndDifferentiate(const Cluster<3> &triplet)
        const {

        Real r01 = triplet.between(0, 1).magnitude;
        Real r02 = triplet.between(0, 2).magnitude;
        Real r12 = triplet.between(1, 2).magnitude;

        auto [f_cut_01, df_cut_01] = cutoff.evaluateAndDifferentiate(r01);
        auto [f_cut_02, df_cut_02] = cutoff.evaluateAndDifferentiate(r02);

        return DescriptorAndDerivatives<4, 3>{
            .value = {
                r01 + r02,
                (r01 - r02) * (r01 - r02),
                r12,
                f_cut_01 * f_cut_02
            },
            .derivatives = std::array{
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
}

#endif