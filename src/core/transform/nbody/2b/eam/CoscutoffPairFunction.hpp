#ifndef JGAP_COSCUTOFFPAIRFUNCTION_HPP
#define JGAP_COSCUTOFFPAIRFUNCTION_HPP

#include <cmath>
#include <tuple>
#include "EamPairFunction.hpp"
#include "core/Real.hpp"

namespace jgap {
    class CoscutoffPairFunction final : public EamPairFunction {
    public:

        CoscutoffPairFunction(const Real cutoff, const Real r_min, const Real prefactor = 1.0)
            : EamPairFunction(cutoff, prefactor), r_min(r_min) {
            interval_inverse = 1.0 / (cutoff - r_min);
        }

        Descriptor<1> evaluate(const Cluster2& pair) const override {
            Real distance = pair.r01;
            if (distance >= cutoff) return {{0.0}};
            if (distance <= r_min) return {{prefactor}};

            const Real chi = (distance - r_min) * interval_inverse;
            return {{prefactor * 0.5 * (1.0 + std::cos(M_PI * chi))}};
        }

        TwoBodyDescriptor<1> evaluateAndDifferentiate(const Cluster2& pair) const override {
            Real distance = pair.r01;
            if (distance >= cutoff) return {{{0.0}}, {}};
            if (distance <= r_min) return {{{prefactor}}, {}};

            const Real chi = (distance - r_min) * interval_inverse;
            const Real dchi_dr = interval_inverse;

            Real val = prefactor * 0.5 * (1.0 + std::cos(M_PI * chi));
            Real deriv = -prefactor * dchi_dr * 0.5 * M_PI * std::sin(M_PI * chi);

            return {{std::array{val}}, {std::array{deriv}}};
        }

        CoscutoffPairFunction* clone() const override {
            return new CoscutoffPairFunction(*this);
        }

        Real getRMin() const {
            return r_min;
        }

    private:
        Real r_min;
        Real interval_inverse;
    };
}

#endif