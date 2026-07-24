#ifndef FSGENPAIRFUNCTION_HPP
#define FSGENPAIRFUNCTION_HPP

#include <cmath>
#include <tuple>
#include "EamPairFunction.hpp"
#include "jgap/core/Real.hpp"

namespace jgap {
    class FSGenPairFunction final : public EamPairFunction {
    public:
        FSGenPairFunction(const Real cutoff, const Real degree, const Real prefactor = 1.0) :
            EamPairFunction(cutoff, prefactor), degree(degree) {
            cutoff_inverse = 1.0 / cutoff;
        }

        Descriptor<1> evaluate(const Cluster2 &pair) const override {
            Real distance = pair.r01;
            if (distance >= cutoff) return {{0.0}};
            return {{prefactor * std::pow(1.0 - distance * cutoff_inverse, degree)}};
        }

        TwoBodyDescriptor<1> evaluateAndDifferentiate(const Cluster2 &pair) const override {
            Real distance = pair.r01;
            if (distance >= cutoff) return {{{0.0}}, {}};

            Real val = prefactor * std::pow(1.0 - distance * cutoff_inverse, degree);
            Real deriv = -prefactor * std::pow(1.0 - distance * cutoff_inverse, degree - 1.0) * degree * cutoff_inverse;

            return {{{val}}, {std::array{deriv}}};
        }

        FSGenPairFunction *clone() const override { return new FSGenPairFunction(*this); }

        Real getDegree() const { return degree; }

    private:
        Real cutoff_inverse;
        Real degree;
    };
}

#endif
