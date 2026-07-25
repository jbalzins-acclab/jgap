#ifndef POLYCUTOFFPAIRFUNCTION_HPP
#define POLYCUTOFFPAIRFUNCTION_HPP

#include <cmath>
#include <tuple>
#include "EamPairFunction.hpp"
#include "jgap/core/Real.hpp"

namespace jgap {
    class PolycutoffPairFunction final : public EamPairFunction {
    public:
        PolycutoffPairFunction(const Real cutoff, const Real r_min, const Real prefactor = 1.0_r) :
            EamPairFunction(cutoff, prefactor), r_min(r_min) {
            interval_inverse = 1.0_r / (cutoff - r_min);
        }

        Descriptor<1> evaluate(const Cluster2& pair) const override {
            Real distance = pair.r01;
            if (distance >= cutoff) return {{0.0_r}};
            if (distance <= r_min) return {{prefactor}};

            const Real chi = (distance - r_min) * interval_inverse;
            return {{prefactor * (1.0_r - chi * chi * chi * (6.0_r * chi * chi - 15.0_r * chi + 10.0_r))}};
        }

        TwoBodyDescriptor<1> evaluateAndDifferentiate(const Cluster2& pair) const override {
            Real distance = pair.r01;
            if (distance >= cutoff) return {.value = {0.0_r}, .derivatives = {0.0_r}};
            if (distance <= r_min) return {.value = {prefactor}, .derivatives = {0.0_r}};

            const Real chi = (distance - r_min) * interval_inverse;
            const Real dchi_dr = interval_inverse;

            Real val = prefactor * (1.0_r - chi * chi * chi * (6.0_r * chi * chi - 15.0_r * chi + 10.0_r));
            Real deriv = prefactor * (dchi_dr * chi * chi * (-30.0_r * chi * chi + 60.0_r * chi - 30.0_r));

            return {.value = {val}, .derivatives = {deriv}};
        }

        PolycutoffPairFunction* clone() const override { return new PolycutoffPairFunction(*this); }

        Real getRMin() const { return r_min; }

    private:
        Real r_min;
        Real interval_inverse;
    };
}

#endif
