#ifndef POLYCUTOFFPAIRFUNCTION_HPP
#define POLYCUTOFFPAIRFUNCTION_HPP

#include <cmath>
#include <tuple>
#include "EamPairFunction.hpp"

namespace jgap {
    class PolycutoffPairFunction final : public EamPairFunction {
    public:
        PolycutoffPairFunction(const double cutoff, const double r_min, const double prefactor = 1.0) :
            EamPairFunction(cutoff, prefactor), r_min(r_min) {
            interval_inverse = 1.0 / (cutoff - r_min);
        }

        Descriptor<1> evaluate(const Cluster2& pair) const override { return TwoBodyTransformation<1>::evaluate(pair); }

        TwoBodyDescriptor<1> evaluateAndDifferentiate(const Cluster2& pair) const override final {
            double distance = pair.separation01.magnitude;
            const auto& dir = pair.separation01.direction;
            if (distance >= cutoff) return {.value = {0.0}, .grad_r1 = {Vector3{}}};
            if (distance <= r_min) return {.value = {prefactor}, .grad_r1 = {Vector3{}}};

            const double chi = (distance - r_min) * interval_inverse;
            const double dchi_dr = interval_inverse;

            double val = prefactor * (1.0 - chi * chi * chi * (6.0 * chi * chi - 15.0 * chi + 10.0));
            double deriv = prefactor * (dchi_dr * chi * chi * (-30.0 * chi * chi + 60.0 * chi - 30.0));

            return {.value = {val}, .grad_r1 = {deriv * dir}};
        }

        PolycutoffPairFunction* clone() const override { return new PolycutoffPairFunction(*this); }
        bool isRotationallyInvariant() const override { return true; }

        double getRMin() const { return r_min; }

    private:
        double r_min;
        double interval_inverse;
    };
}

#endif
