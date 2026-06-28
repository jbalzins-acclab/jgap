#ifndef JGAP_TWOBODYTRANSFORMER_HPP
#define JGAP_TWOBODYTRANSFORMER_HPP

#include "core/cutoff/CutoffFunction.hpp"
#include "../NBodyTransformation.hpp"
#include <memory>

namespace jgap {

    class TwoBodyTransformation final : public NBodyTransformation<2, 2> {
    public:

        TwoBodyTransformation(const ValuePtr<CutoffFunction>& cutoff) : cutoff(cutoff) {
        }

        Descriptor<2> evaluate(const Cluster<2>& pair) const override {
            Real r = pair.separationBetween(0, 1);
            Real f_cut = cutoff->evaluate(r);
            return {{r, f_cut}};
        }

        NBodyDescriptor<2, 2> evaluateAndDifferentiate(const Cluster<2>& pair) const override {
            Real r = pair.separationBetween(0, 1);
            auto [f_cut, df_cut] = cutoff->evaluateAndDifferentiate(r);
            return {{{r, f_cut}}, {std::array{1.0, df_cut}}};
        }

        Cutoffs getCutoffs() const override {
            return Cutoffs{{2, cutoff->getCutoff()}};
        }

        ValuePtr<CutoffFunction> getCutoff() const {
            return cutoff;
        }

        TwoBodyTransformation* clone() const override {
            return new TwoBodyTransformation(*this);
        }

        Real symmetryFactor() const override {
            return 2.0;
        }

    private:
        ValuePtr<CutoffFunction> cutoff;
    };
}

#endif