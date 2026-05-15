#ifndef JGAP_ANGLE3BTRANSFORMER_HPP
#define JGAP_ANGLE3BTRANSFORMER_HPP

#include <utility>

#include "core/cutoff/CutoffFunction.hpp"
#include "../Transformer.hpp"
#include "core/transform/KnownDependenciesTransformer.hpp"

namespace jgap {

    template<typename TCutoff>
    requires std::derived_from<TCutoff, CutoffFunction>
    class Angle3bTransformer : public KnownDependenciesTransformer<3+1, 3> {
    public:
        ~Angle3bTransformer() override = default;

        Angle3bTransformer() = default;
        Angle3bTransformer(TCutoff cutoff) : cutoff(cutoff) { }

        std::array<Real, 4> value(const Separations<3> &separations) const override;
        ValueAndDerivatives valueAndDerivatives(const Separations<3> &separations) const override;

        Real getCutoff() const override { return cutoff.getCutoff(); }

    private:
        TCutoff cutoff;
    };

    template<typename TCutoff> requires std::derived_from<TCutoff, CutoffFunction>
    std::array<Real, 4> Angle3bTransformer<TCutoff>::value(const Separations<3> &separations) const {

        Real r01 = separations.between(0, 1).magnitude;
        Real r02 = separations.between(0, 2).magnitude;
        Real r12 = separations.between(1, 2).magnitude;

        Real f_cut_01 = cutoff.evaluate(r01);
        Real f_cut_02 = cutoff.evaluate(r02);

        return {
            r01 + r02,
            (r01 - r02) * (r01 - r02),
            r12,
            f_cut_01 * f_cut_02
        };
    }

    template<typename TCutoff> requires std::derived_from<TCutoff, CutoffFunction>
    KnownDependenciesTransformer<4, 3>::ValueAndDerivatives Angle3bTransformer<TCutoff>::valueAndDerivatives(
        const Separations<3> &separations) const {

        Real r01 = separations.between(0, 1).magnitude;
        Real r02 = separations.between(0, 2).magnitude;
        Real r12 = separations.between(1, 2).magnitude;

        Real f_cut_01 = cutoff.evaluate(r01);
        Real f_cut_02 = cutoff.evaluate(r02);

        Real df_cut_01 = cutoff.differentiate(r01);
        Real df_cut_02 = cutoff.differentiate(r02);

        return {
            .value = std::array{
                r01 + r02,
                (r01 - r02) * (r01 - r02),
                r12,
                f_cut_01 * f_cut_02
            },
            .derivatives = std::array{
                std::array{ // wrt r_01
                    static_cast<Real>(1.0),
                    static_cast<Real>(2.0) * (r01 - r02),
                    static_cast<Real>(0.0),
                    df_cut_01 * f_cut_02
                },
                std::array{ // wrt r_02
                    static_cast<Real>(1.0),
                    static_cast<Real>(2.0) * (r02 - r01),
                    static_cast<Real>(0.0),
                    df_cut_02 * f_cut_01
                },
                std::array{ // wrt r_12
                    static_cast<Real>(0.0),
                    static_cast<Real>(0.0),
                    static_cast<Real>(1.0),
                    static_cast<Real>(0.0)
                }
            }
        };
    }
}

#endif