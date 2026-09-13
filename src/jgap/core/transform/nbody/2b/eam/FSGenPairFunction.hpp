#ifndef FSGENPAIRFUNCTION_HPP
#define FSGENPAIRFUNCTION_HPP

#include <cmath>
#include <tuple>
#include "EamPairFunction.hpp"

namespace jgap {
    class FSGenPairFunction final : public EamPairFunction {
    public:
        FSGenPairFunction(const double cutoff, const double degree, const double prefactor = 1.0) :
            EamPairFunction(cutoff, prefactor), degree(degree) {
            cutoff_inverse = 1.0 / cutoff;
        }

        Descriptor<1> evaluate(const Cluster2& pair) const override { return TwoBodyTransformation<1>::evaluate(pair); }

        TwoBodyDescriptor<1> evaluateAndDifferentiate(const Cluster2& pair) const override final {
            double distance = pair.separation01.magnitude;
            const auto& dir = pair.separation01.direction;
            if (distance >= cutoff) return { .value = { 0.0 }, .grad_r1 = { Vector3{} } };

            double val = prefactor * std::pow(1.0 - distance * cutoff_inverse, degree);
            double deriv =
                -prefactor * std::pow(1.0 - distance * cutoff_inverse, degree - 1.0) * degree * cutoff_inverse;

            return { .value = { val }, .grad_r1 = { deriv * dir } };
        }

        FSGenPairFunction* clone() const override { return new FSGenPairFunction(*this); }
        bool isRotationallyInvariant() const override { return true; }

        double getDegree() const { return degree; }

    private:
        double cutoff_inverse;
        double degree;
    };
}

#endif
