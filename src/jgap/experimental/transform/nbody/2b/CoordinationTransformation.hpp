#ifndef JGAP_COORDINATIONSTRANSFORMATION_HPP
#define JGAP_COORDINATIONSTRANSFORMATION_HPP

#include "jgap/core/transform/nbody/2b/TwoBodyTransformation.hpp"
#include "jgap/experimental/cutoff/WendlandFunction.hpp"
#include <array>
#include <algorithm>
#include <utility>

namespace jgap {

    template<size_t Dim>
    class CoordinationTransformation final : public TwoBodyTransformation<Dim> {
    public:
        CoordinationTransformation(const std::array<std::pair<Real, Real>, Dim>& ranges) : ranges(ranges) {
            max_cutoff = 0.0_r;
            for (size_t i = 0; i < Dim; ++i) {
                wendlands[i] = WendlandFunction(ranges[i].first, ranges[i].second);
                if (ranges[i].second > max_cutoff) {
                    max_cutoff = ranges[i].second;
                }
            }
        }

        Descriptor<Dim> evaluate(const Cluster2& pair) const override { 
            return TwoBodyTransformation<Dim>::evaluate(pair); 
        }

        TwoBodyDescriptor<Dim> evaluateAndDifferentiate(const Cluster2& pair) const override final {
            Real r = pair.separation01.magnitude;
            TwoBodyDescriptor<Dim> desc;
            for (size_t i = 0; i < Dim; ++i) {
                auto [val, deriv] = wendlands[i].evaluateAndDifferentiate(r);
                desc.value[i] = val;
                desc.derivatives[i] = deriv;
            }
            return desc;
        }

        Cutoffs getCutoffs() const override { return Cutoffs{{2, max_cutoff}}; }

        const std::array<std::pair<Real, Real>, Dim>& getRanges() const { return ranges; }

        CoordinationTransformation* clone() const override { return new CoordinationTransformation(*this); }

    private:
        std::array<std::pair<Real, Real>, Dim> ranges;
        std::array<WendlandFunction, Dim> wendlands;
        Real max_cutoff;
    };
}

#endif // JGAP_COORDINATIONSTRANSFORMATION_HPP
