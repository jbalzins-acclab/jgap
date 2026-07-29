#ifndef JGAP_PERRIOTPOLYNOMIALCUTOFF_HPP
#define JGAP_PERRIOTPOLYNOMIALCUTOFF_HPP

#include <tuple>
#include "CutoffFunction.hpp"
#include "jgap/core/Real.hpp"

namespace jgap {
    class PerriotPolynomialCutoff final : public CutoffFunction {
    public:
        PerriotPolynomialCutoff(Real cutoff, Real cutoff_transition_width) :
            cutoff(cutoff),
            r_min(cutoff - cutoff_transition_width),
            cutoff_width_inverse(1.0_r / cutoff_transition_width) {}

        Real getCutoff() const override { return cutoff; }
        Real getCutoffTransitionWidth() const { return cutoff - r_min; }

        Real evaluate(Real r) const override {
            return CutoffFunction::evaluate(r);
        }

        std::tuple<Real, Real> evaluateAndDifferentiate(Real r) const override {
            if (r <= r_min) return {1.0_r, 0.0_r};
            if (r >= cutoff) [[unlikely]]
                return {0.0_r, 0.0_r};

            const Real chi = (r - r_min) * cutoff_width_inverse;
            const Real chi_sq = chi * chi;
            const Real chi_cube = chi_sq * chi;
            Real val = 1.0_r - chi_cube * (10.0_r - 15.0_r * chi + 6.0_r * chi_sq);
            Real deriv = -30.0_r * chi_sq * (1.0_r - 2.0_r * chi + chi_sq) * cutoff_width_inverse;

            return {val, deriv};
        }

        PerriotPolynomialCutoff* clone() const override { return new PerriotPolynomialCutoff(*this); }

    private:
        const Real cutoff;
        const Real r_min;
        const Real cutoff_width_inverse;
    };
}

#endif
