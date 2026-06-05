#ifndef JGAP_PERRIOTPOLYNOMIALCUTOFF_HPP
#define JGAP_PERRIOTPOLYNOMIALCUTOFF_HPP

#include <tuple>
#include "CutoffFunction.hpp"
#include "core/Real.hpp"

namespace jgap {
    class PerriotPolynomialCutoff final : public CutoffFunction {
    public:
        PerriotPolynomialCutoff(Real r_min, Real cutoff)
            : cutoff(cutoff),
              r_min(r_min),
              cutoff_width_inverse(1.0 / (cutoff - r_min)) {}

        Real getCutoff() const override { return cutoff; }

        Real evaluate(Real r) const override {
            if (r <= r_min) return 1.0;
            if (r >= cutoff) return 0.0;

            const Real chi = (r - r_min) * cutoff_width_inverse;
            return 1.0 - chi * chi * chi * (10.0 - 15.0 * chi + 6.0 * chi * chi);
        }

        std::tuple<Real, Real> evaluateAndDifferentiate(Real r) const override {
            if (r <= r_min) return {1.0, 0.0};
            if (r >= cutoff) return {0.0, 0.0};

            const Real chi = (r - r_min) * cutoff_width_inverse;
            const Real chi_sq = chi * chi;
            const Real chi_cube = chi * chi;
            Real val = 1.0 - chi_cube * (10.0 - 15.0 * chi + 6.0 * chi_sq);
            Real deriv = -30.0 * chi_sq * (1.0 - 2.0 * chi + chi_sq) * cutoff_width_inverse;

            return {val, deriv};
        }

    private:
        const Real cutoff;
        const Real r_min;
        const Real cutoff_width_inverse;
    };
}

#endif