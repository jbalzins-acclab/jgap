#ifndef JGAP_PAIRDISTANCETRANSFORMATION_HPP
#define JGAP_PAIRDISTANCETRANSFORMATION_HPP

#include "TwoBodyTransformation.hpp"
#include "jgap/core/cutoff/CutoffFunction.hpp"

namespace jgap {

    class PairDistanceTransformation final : public TwoBodyTransformation<2> {
    public:
        PairDistanceTransformation(const ValuePtr<CutoffFunction> &cutoff) : cutoff(cutoff) {}

        Descriptor<2> evaluate(const Cluster2 &pair) const override {
            Real r = pair.r01;
            Real f_cut = cutoff->evaluate(r);
            return {{r, f_cut}};
        }

        TwoBodyDescriptor<2> evaluateAndDifferentiate(const Cluster2 &pair) const override {
            Real r = pair.r01;
            auto [f_cut, df_cut] = cutoff->evaluateAndDifferentiate(r);
            return {{{r, f_cut}}, {std::array{1.0, df_cut}}};
        }

        Cutoffs getCutoffs() const override { return Cutoffs{{2, cutoff->getCutoff()}}; }

        ValuePtr<CutoffFunction> getCutoffFunction() const { return cutoff; }

        PairDistanceTransformation *clone() const override { return new PairDistanceTransformation(*this); }

    private:
        ValuePtr<CutoffFunction> cutoff;
    };
}

#endif
