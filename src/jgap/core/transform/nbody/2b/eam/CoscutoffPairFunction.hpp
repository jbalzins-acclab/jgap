#ifndef JGAP_COSCUTOFFPAIRFUNCTION_HPP
#define JGAP_COSCUTOFFPAIRFUNCTION_HPP

#include <cmath>
#include <tuple>
#include "EamPairFunction.hpp"

namespace jgap {
    class CoscutoffPairFunction final : public EamPairFunction {
    public:
        CoscutoffPairFunction(const double cutoff, const double r_min, const double prefactor = 1.0) :
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

            double val = prefactor * 0.5 * (1.0 + std::cos(static_cast<double>(M_PI) * chi));
            double deriv =
                -prefactor * dchi_dr * 0.5 * static_cast<double>(M_PI) * std::sin(static_cast<double>(M_PI) * chi);

            return {.value = {val}, .grad_r1 = {deriv * dir}};
        }

        CoscutoffPairFunction* clone() const override { return new CoscutoffPairFunction(*this); }
        bool isRotationallyInvariant() const override { return true; }

        double getRMin() const { return r_min; }

    private:
        double r_min;
        double interval_inverse;
    };
}

#endif
