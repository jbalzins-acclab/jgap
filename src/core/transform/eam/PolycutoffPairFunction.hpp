#ifndef POLYCUTOFFPAIRFUNCTION_HPP
#define POLYCUTOFFPAIRFUNCTION_HPP

#include <cmath>
#include <tuple>
#include "EamPairFunction.hpp"
#include "core/Real.hpp"

namespace jgap {
    class PolycutoffPairFunction final : public EamPairFunction {
    public:
        PolycutoffPairFunction(const Real cutoff, const Real r_min, const Real prefactor = 1.0)
            : EamPairFunction(cutoff, prefactor), r_min(r_min)
        {
            interval_inverse = 1.0 / (cutoff - r_min);
        }

        Descriptor<1> evaluate(const Cluster<2>& pair) const override {
            Real distance = pair.between(0, 1).magnitude;
            if (distance >= cutoff) return {{0.0}};
            if (distance <= r_min) return {{prefactor}};

            const Real chi = (distance - r_min) * interval_inverse;
            return {{prefactor * (1.0 - chi * chi * chi * (6.0 * chi * chi - 15.0 * chi + 10.0))}};
        }

        NBodyDescriptor<1, 2> evaluateAndDifferentiate(const Cluster<2>& pair) const override {
            Real distance = pair.between(0, 1).magnitude;
            if (distance >= cutoff) return {{{0.0}}, {}};
            if (distance <= r_min) return {{{prefactor}}, {}};

            const Real chi = (distance - r_min) * interval_inverse;
            const Real dchi_dr = interval_inverse;

            Real val = prefactor * (1.0 - chi * chi * chi * (6.0 * chi * chi - 15.0 * chi + 10.0));
            Real deriv = prefactor * (dchi_dr * chi * chi * ( -30.0 * chi * chi + 60.0 * chi - 30.0));

            return {{{val}}, {std::array{deriv}}};
        }

    private:
        Real r_min;
        Real interval_inverse;
    };
}

#endif