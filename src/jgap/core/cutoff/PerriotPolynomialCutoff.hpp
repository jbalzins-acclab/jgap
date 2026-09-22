#ifndef JGAP_PERRIOTPOLYNOMIALCUTOFF_HPP
#define JGAP_PERRIOTPOLYNOMIALCUTOFF_HPP

#include <tuple>
#include "CutoffFunction.hpp"

namespace jgap {
    class PerriotPolynomialCutoff final : public CutoffFunction {
    public:
        PerriotPolynomialCutoff(double cutoff, double cutoff_transition_width) :
            cutoff(cutoff),
            r_min(cutoff - cutoff_transition_width),
            cutoff_width_inverse(1.0 / cutoff_transition_width) {}

        double getCutoff() const override { return cutoff; }
        double getCutoffTransitionWidth() const { return cutoff - r_min; }

        double evaluate(double r) const override {
            return CutoffFunction::evaluate(r);
        }

        std::tuple<double, double> evaluateAndDifferentiate(double r) const override {
            if (r <= r_min) return {1.0, 0.0};
            if (r >= cutoff) [[unlikely]]
                return {0.0, 0.0};

            const double chi = (r - r_min) * cutoff_width_inverse;
            const double chi_sq = chi * chi;
            const double chi_cube = chi_sq * chi;
            double val = 1.0 - chi_cube * (10.0 - 15.0 * chi + 6.0 * chi_sq);
            double deriv = -30.0 * chi_sq * (1.0 - 2.0 * chi + chi_sq) * cutoff_width_inverse;

            return {val, deriv};
        }

        PerriotPolynomialCutoff* clone() const override { return new PerriotPolynomialCutoff(*this); }

    private:
        const double cutoff;
        const double r_min;
        const double cutoff_width_inverse;
    };
}

#endif
