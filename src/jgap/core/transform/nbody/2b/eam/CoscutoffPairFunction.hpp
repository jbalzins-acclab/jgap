#ifndef JGAP_COSCUTOFFPAIRFUNCTION_HPP
#define JGAP_COSCUTOFFPAIRFUNCTION_HPP

#include <cmath>
#include <tuple>
#include "EamPairFunction.hpp"
#include "jgap/core/Real.hpp"

namespace jgap {
    class CoscutoffPairFunction final : public EamPairFunction {
    public:
        CoscutoffPairFunction(const Real cutoff, const Real r_min, const Real prefactor = 1.0_r) :
            EamPairFunction(cutoff, prefactor), r_min(r_min) {
            interval_inverse = 1.0_r / (cutoff - r_min);
        }

        Descriptor<1> evaluate(const Cluster2& pair) const override { return TwoBodyTransformation<1>::evaluate(pair); }

        TwoBodyDescriptor<1> evaluateAndDifferentiate(const Cluster2& pair) const override final {
            Real distance = pair.separation01.magnitude;
            const auto& dir = pair.separation01.direction;
            if (distance >= cutoff) return {.value = {0.0_r}, .grad_r1 = {Vector3{}}};
            if (distance <= r_min) return {.value = {prefactor}, .grad_r1 = {Vector3{}}};

            const Real chi = (distance - r_min) * interval_inverse;
            const Real dchi_dr = interval_inverse;

            Real val = prefactor * 0.5_r * (1.0_r + std::cos(static_cast<Real>(M_PI) * chi));
            Real deriv =
                -prefactor * dchi_dr * 0.5_r * static_cast<Real>(M_PI) * std::sin(static_cast<Real>(M_PI) * chi);

            return {.value = {val}, .grad_r1 = {deriv * dir}};
        }

        CoscutoffPairFunction* clone() const override { return new CoscutoffPairFunction(*this); }
        bool isRotationallyInvariant() const override { return true; }

        Real getRMin() const { return r_min; }

    private:
        Real r_min;
        Real interval_inverse;
    };
}

#endif
