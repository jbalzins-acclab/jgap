#ifndef FSGENPAIRFUNCTION_HPP
#define FSGENPAIRFUNCTION_HPP

#include <cmath>
#include <tuple>
#include "EamPairFunction.hpp"
#include "jgap/core/Real.hpp"

namespace jgap {
    class FSGenPairFunction final : public EamPairFunction {
    public:
        FSGenPairFunction(const Real cutoff, const Real degree, const Real prefactor = 1.0_r) :
            EamPairFunction(cutoff, prefactor), degree(degree) {
            cutoff_inverse = 1.0_r / cutoff;
        }

        Descriptor<1> evaluate(const Cluster2& pair) const override {
            return TwoBodyTransformation<1>::evaluate(pair);
        }

        TwoBodyDescriptor<1> evaluateAndDifferentiate(const Cluster2& pair) const override final {
            Real distance = pair.separation01.magnitude;
            if (distance >= cutoff) return {.value = {0.0_r}, .derivatives = {0.0_r}};

            Real val = prefactor * std::pow(1.0_r - distance * cutoff_inverse, degree);
            Real deriv =
                -prefactor * std::pow(1.0_r - distance * cutoff_inverse, degree - 1.0_r) * degree * cutoff_inverse;

            return {.value = {val}, .derivatives = {deriv}};
        }

        FSGenPairFunction* clone() const override { return new FSGenPairFunction(*this); }

        Real getDegree() const { return degree; }

    private:
        Real cutoff_inverse;
        Real degree;
    };
}

#endif
