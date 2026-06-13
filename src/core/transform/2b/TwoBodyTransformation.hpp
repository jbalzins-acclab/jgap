#ifndef JGAP_TWOBODYTRANSFORMER_HPP
#define JGAP_TWOBODYTRANSFORMER_HPP

#include "core/cutoff/CutoffFunction.hpp"
#include "../ClusterTransformation.hpp"
#include <memory>

namespace jgap {

    class TwoBodyTransformation final : public ClusterTransformation<2, 2> {
    public:

        TwoBodyTransformation(const ValuePtr<CutoffFunction>& cutoff) : cutoff(cutoff) {
        }

        Descriptor<2> evaluate(const Cluster<2>& pair) const override {
            Real r = pair.between(0, 1);
            Real f_cut = cutoff->evaluate(r);
            return {{r, f_cut}};
        }

        NBodyDescriptor<2, 2> evaluateAndDifferentiate(const Cluster<2>& pair) const override {
            Real r = pair.between(0, 1);
            auto [f_cut, df_cut] = cutoff->evaluateAndDifferentiate(r);
            return {{{r, f_cut}}, {std::array{1.0, df_cut}}};
        }

        Cutoffs getCutoffs() const override {
            return Cutoffs{{2, cutoff->getCutoff()}};
        }

        std::unique_ptr<ClusterTransformation> clone() const override {
            return std::make_unique<TwoBodyTransformation>(*this);
        }

    private:
        ValuePtr<CutoffFunction> cutoff;
    };
}

#endif