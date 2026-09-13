#ifndef JGAP_COSCUTOFF_HPP
#define JGAP_COSCUTOFF_HPP

#include <cmath>
#include <tuple>
#include "CutoffFunction.hpp"
#include "jgap/utils/Utils.hpp"

namespace jgap {
    class CosCutoff final : public CutoffFunction {
    public:
        CosCutoff(double cutoff, double cutoff_transition_width) :
            cutoff(cutoff),
            r_min(cutoff - cutoff_transition_width),
            pi_over_w(static_cast<double>(M_PI) / cutoff_transition_width),
            deriv_coeff(-0.5 * pi_over_w) {}

        double getCutoff() const override { return cutoff; }
        double getCutoffTransitionWidth() const { return cutoff - r_min; }

        double evaluate(double r) const override { return CutoffFunction::evaluate(r); }

        std::tuple<double, double> evaluateAndDifferentiate(double r) const override final {
            if (r <= r_min) return {1.0, 0.0};
            if (r >= cutoff) [[unlikely]]
                return {0.0, 0.0};

            const double phase = (r - r_min) * pi_over_w;
            double s, c;
            utils::sincos(phase, &s, &c);

            double val = 0.5 * (c + 1.0);
            double deriv = deriv_coeff * s;

            return {val, deriv};
        }

        CosCutoff* clone() const override { return new CosCutoff(*this); }

    private:
        const double cutoff;
        const double r_min;
        const double pi_over_w;
        const double deriv_coeff;
    };
}
#endif
