#ifndef JGAP_SPLINEPAIRTRANSFORMATION_HPP
#define JGAP_SPLINEPAIRTRANSFORMATION_HPP
#include <utility>

#include "../TwoBodyTransformation.hpp"
#include "jgap/core/splines/HermiteCubicSpline.hpp"

namespace jgap {
    class SplinePairTransformation final : public TwoBodyTransformation<1> {
    public:
        SplinePairTransformation(const SplinePairTransformation& other) = default;
        SplinePairTransformation(HermiteCubicSpline spline) : spline(std::move(spline)) {}

        Descriptor<1> evaluate(const Cluster2& cluster) const override {
            return TwoBodyTransformation<1>::evaluate(cluster);
        }

        TwoBodyDescriptor<1> evaluateAndDifferentiate(const Cluster2& cluster) const override final {
            Real r01 = cluster.separation01.magnitude;
            const auto& dir = cluster.separation01.direction;

            auto [val, derivative] = spline.interpolate({ r01 });

            return { .value = { val }, .grad_r1 = { derivative[0] * dir } };
        }

        Cutoffs getCutoffs() const override { return { { 2u, spline.getCutoff()[0] } }; }
        bool isRotationallyInvariant() const override { return true; }

        SplinePairTransformation* clone() const override { return new SplinePairTransformation(*this); }

        const auto& getSpline() const { return spline; }

    private:
        HermiteCubicSpline spline;
    };
}

#endif
