#ifndef JGAP_WENDLANDFUNCTION_HPP
#define JGAP_WENDLANDFUNCTION_HPP

#include "jgap/core/cutoff/CutoffFunction.hpp"
#include "../../core/io/log/CurrentLogger.hpp"

namespace jgap {

    class WendlandFunction final : public CutoffFunction {
    public:
        WendlandFunction() : r_min(0.0), r_max(1.0), inv_range(1.0) {}

        WendlandFunction(double r_min, double r_max) : r_min(r_min), r_max(r_max) {
            if (r_max <= r_min) {
                JGAP_LOG_AND_THROW("WendlandFunction requires r_max > r_min");
            }
            inv_range = 1.0 / (r_max - r_min);
        }

        std::tuple<double, double> evaluateAndDifferentiate(double r) const override {
            if (r <= r_min || r >= r_max) {
                return {0.0, 0.0};
            }

            double x = (r - r_min) * inv_range;
            double omx = 1.0 - x;
            double omx2 = omx * omx;
            double omx3 = omx2 * omx;
            double omx4 = omx2 * omx2;

            double val = omx4 * (1.0 + 4.0 * x);
            double grad_x = -20.0 * x * omx3;
            double grad_r = grad_x * inv_range;

            return {val, grad_r};
        }

        double getCutoff() const override { return r_max; }

        double getRMin() const { return r_min; }
        double getRMax() const { return r_max; }

        WendlandFunction* clone() const override { return new WendlandFunction(*this); }

    private:
        double r_min;
        double r_max;
        double inv_range;
    };
}

#endif // JGAP_WENDLANDFUNCTION_HPP
