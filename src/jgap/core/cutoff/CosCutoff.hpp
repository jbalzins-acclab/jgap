#ifndef JGAP_COSCUTOFF_HPP
#define JGAP_COSCUTOFF_HPP

#include <cmath>
#include <tuple>
#include "CutoffFunction.hpp"
#include "jgap/core/Real.hpp"
#include "jgap/utils/Utils.hpp"

namespace jgap {
    class CosCutoff final : public CutoffFunction {
    public:
        CosCutoff(Real cutoff, Real cutoff_transition_width) :
            cutoff(cutoff),
            r_min(cutoff - cutoff_transition_width),
            pi_over_w(static_cast<Real>(M_PI) / cutoff_transition_width),
            deriv_coeff(-0.5_r * pi_over_w) {}

        Real getCutoff() const override { return cutoff; }
        Real getCutoffTransitionWidth() const { return cutoff - r_min; }

        Real evaluate(Real r) const override {
            return CutoffFunction::evaluate(r);
        }

        std::tuple<Real, Real> evaluateAndDifferentiate(Real r) const override final {
            if (r <= r_min) return {1.0_r, 0.0_r};
            if (r >= cutoff) [[unlikely]]
                return {0.0_r, 0.0_r};

            const Real phase = (r - r_min) * pi_over_w;
            Real s, c;
            jgap::sincos(phase, &s, &c);

            Real val = 0.5_r * (c + 1.0_r);
            Real deriv = deriv_coeff * s;

            return {val, deriv};
        }

        CosCutoff* clone() const override { return new CosCutoff(*this); }

    private:
        const Real cutoff;
        const Real r_min;
        const Real pi_over_w;
        const Real deriv_coeff;
    };
}
#endif
