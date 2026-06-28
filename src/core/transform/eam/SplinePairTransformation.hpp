#ifndef JGAP_SPLINEPAIRTRANSFORMATION_HPP
#define JGAP_SPLINEPAIRTRANSFORMATION_HPP
#include <utility>

#include "core/splines/Spline.hpp"
#include "core/transform/NBodyTransformation.hpp"

namespace jgap {
    class SplinePairTransformation : public NBodyTransformation<1, 2> {
    public:

        SplinePairTransformation(const SplinePairTransformation& other) = default;
        SplinePairTransformation(HermiteCubicSpline spline) : spline(std::move(spline)) {}

        Descriptor<1> evaluate(const Cluster<2> &cluster) const override {
            Real r01 = cluster.separationBetween(0, 1);

            return {
                .value = {
                    spline.interpolate({r01}).value
                }
            };
        }

        NBodyDescriptor<1, 2> evaluateAndDifferentiate(const Cluster<2> &cluster) const override {
            Real r01 = cluster.separationBetween(0, 1);

            auto [val, derivative] = spline.interpolate({r01});

            return {
                .value = {
                    val
                },
                .derivatives = {
                    derivative
                }
            };
        }

        Cutoffs getCutoffs() const override {
            return {{2u, spline.getCutoff()[0]}};
        }

        SplinePairTransformation* clone() const override {
            return new SplinePairTransformation(*this);
        }

        const auto& getSpline() const {
            return spline;
        }

    private:
        HermiteCubicSpline spline;
    };
}

#endif
