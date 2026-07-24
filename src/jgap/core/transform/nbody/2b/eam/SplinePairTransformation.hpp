#ifndef JGAP_SPLINEPAIRTRANSFORMATION_HPP
#define JGAP_SPLINEPAIRTRANSFORMATION_HPP
#include <utility>

#include "../TwoBodyTransformation.hpp"
#include "jgap/core/splines/HermiteCubicSpline.hpp"

namespace jgap {
    class SplinePairTransformation final : public TwoBodyTransformation<1> {
    public:
        SplinePairTransformation(const SplinePairTransformation &other) = default;
        SplinePairTransformation(HermiteCubicSpline spline) : spline(std::move(spline)) {}

        Descriptor<1> evaluate(const Cluster2 &cluster) const override {
            Real r01 = cluster.r01;

            return {spline.interpolate({r01}).value};
        }

        TwoBodyDescriptor<1> evaluateAndDifferentiate(const Cluster2 &cluster) const override {
            Real r01 = cluster.r01;

            auto [val, derivative] = spline.interpolate({r01});

            return {.value = {val}, .derivatives = {derivative}};
        }

        Cutoffs getCutoffs() const override { return {{2u, spline.getCutoff()[0]}}; }

        SplinePairTransformation *clone() const override { return new SplinePairTransformation(*this); }

        const auto &getSpline() const { return spline; }

    private:
        HermiteCubicSpline spline;
    };
}

#endif
