#ifndef JGAP_PAIRDISTANCETRANSFORMATION_HPP
#define JGAP_PAIRDISTANCETRANSFORMATION_HPP

#include "TwoBodyTransformation.hpp"
#include "jgap/core/cutoff/CutoffFunction.hpp"

namespace jgap {

    class PairDistanceTransformation final : public TwoBodyTransformation<2> {
    public:
        PairDistanceTransformation(const ValuePtr<CutoffFunction>& cutoff) : cutoff(cutoff) {}

        Descriptor<2> evaluate(const Cluster2& pair) const override { return TwoBodyTransformation<2>::evaluate(pair); }

        TwoBodyDescriptor<2> evaluateAndDifferentiate(const Cluster2& pair) const override final {
            Real r = pair.separation01.magnitude;
            const auto& dir = pair.separation01.direction;
            auto [f_cut, df_cut] = cutoff->evaluateAndDifferentiate(r);
            return {
                .value = {r, f_cut},
                .grad_r1 = {dir, df_cut * dir},
            };
        }

        Cutoffs getCutoffs() const override { return Cutoffs{{2, cutoff->getCutoff()}}; }
        bool isRotationallyInvariant() const override { return true; }

        ValuePtr<CutoffFunction> getCutoffFunction() const { return cutoff; }

        PairDistanceTransformation* clone() const override { return new PairDistanceTransformation(*this); }

    private:
        ValuePtr<CutoffFunction> cutoff;
    };
}

#endif
