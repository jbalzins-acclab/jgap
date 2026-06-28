#ifndef JGAP_COSCUTOFF_HPP
#define JGAP_COSCUTOFF_HPP

#include <cmath>
#include <tuple>
#include "CutoffFunction.hpp"
#include "core/Real.hpp"
#include "utils/Utils.hpp"

namespace jgap {
    class CosCutoff final : public CutoffFunction {
    public:

        CosCutoff(Real cutoff, Real cutoff_transition_width)
            : cutoff(cutoff),
              r_min(cutoff - cutoff_transition_width),
              pi_over_w(M_PI / cutoff_transition_width),
              deriv_coeff(-0.5 * pi_over_w) {}

        Real getCutoff() const override { return cutoff; }
        Real getCutoffTransitionWidth() const { return cutoff - r_min; }

        Real evaluate(Real r) const override {
            if (r <= r_min) return 1.0;
            if (r >= cutoff) [[unlikely]] return 0.0;
            return 0.5 * (std::cos((r - r_min) * pi_over_w) + 1.0);
        }

        std::tuple<Real, Real> evaluateAndDifferentiate(Real r) const override {
            if (r <= r_min) return {1.0, 0.0};
            if (r >= cutoff) [[unlikely]] return {0.0, 0.0};

            const Real phase = (r - r_min) * pi_over_w;
            Real s, c;
            jgap::sincos(phase, &s, &c);

            Real val = 0.5 * (c + 1.0);
            Real deriv = deriv_coeff * s;

            return {val, deriv};
        }

        CosCutoff* clone() const override {
            return new CosCutoff(*this);
        }

    private:
        const Real cutoff;
        const Real r_min;
        const Real pi_over_w;
        const Real deriv_coeff;
    };
}
#endif