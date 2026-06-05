#ifndef JGAP_TWOBODYTRANSFORMER_HPP
#define JGAP_TWOBODYTRANSFORMER_HPP

#include "core/cutoff/CutoffFunction.hpp"
#include "../ClusterTransformation.hpp"

namespace jgap {

    template<typename TCutoff = CosCutoff>
    requires std::derived_from<TCutoff, CutoffFunction>
    class TwoBodyTransformation : public ClusterTransformation<2, 2> {
    public:
        static constexpr size_t Dim = 2;
        static constexpr size_t ClusterSize = 2;

        TwoBodyTransformation() = default;
        TwoBodyTransformation(TCutoff cutoff) : cutoff(cutoff) { }

        Descriptor<2> evaluate(const Cluster<2>& pair) const override;

        DescriptorAndDerivatives<2, 2> evaluateAndDifferentiate(const Cluster<2>& pair) const override;

        Cutoffs getCutoffs() const override {
            return Cutoffs{ {2, cutoff.getCutoff()} };
        }

    private:
        TCutoff cutoff;
    };

    static_assert(CClusterTransformation<TwoBodyTransformation<>>);

    template<typename TCutoff> requires std::derived_from<TCutoff, CutoffFunction>
    Descriptor<2> TwoBodyTransformation<TCutoff>::evaluate(const Cluster<2>& pair) const {
        Real r = pair.between(0, 1).magnitude;
        Real f_cut = cutoff.evaluate(r);

        return {
            .value = { r, f_cut }
        };
    }

    template<typename TCutoff> requires std::derived_from<TCutoff, CutoffFunction>
    DescriptorAndDerivatives<2, 2> TwoBodyTransformation<TCutoff>::evaluateAndDifferentiate(const Cluster<2>& pair)
        const {

        Real r = pair.between(0, 1).magnitude;
        auto [f_cut, df_cut] = cutoff.evaluateAndDifferentiate(r);

        return {
            .value = { r, f_cut },
            .derivatives = {
                std::array{ 1.0, df_cut }
            }
        };
    }
}

#endif