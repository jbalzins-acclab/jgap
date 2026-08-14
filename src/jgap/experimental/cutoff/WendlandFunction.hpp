#ifndef JGAP_WENDLANDFUNCTION_HPP
#define JGAP_WENDLANDFUNCTION_HPP

#include "jgap/core/cutoff/CutoffFunction.hpp"
#include "../../core/io/log/CurrentLogger.hpp"

namespace jgap {

    class WendlandFunction final : public CutoffFunction {
    public:
        WendlandFunction() : r_min(0.0_r), r_max(1.0_r), inv_range(1.0_r) {}

        WendlandFunction(Real r_min, Real r_max) : r_min(r_min), r_max(r_max) {
            if (r_max <= r_min) {
                JGAP_LOG_AND_THROW("WendlandFunction requires r_max > r_min");
            }
            inv_range = 1.0_r / (r_max - r_min);
        }

        std::tuple<Real, Real> evaluateAndDifferentiate(Real r) const override {
            if (r <= r_min || r >= r_max) {
                return {0.0_r, 0.0_r};
            }

            Real x = (r - r_min) * inv_range;
            Real omx = 1.0_r - x;
            Real omx2 = omx * omx;
            Real omx3 = omx2 * omx;
            Real omx4 = omx2 * omx2;

            Real val = omx4 * (1.0_r + 4.0_r * x);
            Real grad_x = -20.0_r * x * omx3;
            Real grad_r = grad_x * inv_range;

            return {val, grad_r};
        }

        Real getCutoff() const override { return r_max; }

        Real getRMin() const { return r_min; }
        Real getRMax() const { return r_max; }

        WendlandFunction* clone() const override { return new WendlandFunction(*this); }

    private:
        Real r_min;
        Real r_max;
        Real inv_range;
    };
}

#endif // JGAP_WENDLANDFUNCTION_HPP
