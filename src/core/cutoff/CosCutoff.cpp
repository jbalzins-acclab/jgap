#include "core/cutoff/CosCutoff.hpp"

#include <cmath>

#include "utils/Utils.hpp"

namespace jgap {

    CosCutoff::CosCutoff(const double cutoff, const double cutoffTransitionWidth) {
        _cutoff = cutoff;
        _cutoffTransitionWidth = cutoffTransitionWidth;
        _cutoffTransitionWidthInverse = 1.0 / _cutoffTransitionWidth;
    }

    shared_ptr<CosCutoff> CosCutoff::fromJson(nlohmann::json params) {

        double cutoff = require(params, "cutoff");
        double cutoffTransitionWidth;
        if (params.contains("cutoff_transition_width")) {
            cutoffTransitionWidth = params["cutoff_transition_width"];
        } else {
            cutoffTransitionWidth = cutoff - require(params, "r_min").get<double>();
        }

        return make_shared<CosCutoff>(cutoff, cutoffTransitionWidthInverse);
    }

    nlohmann::json CosCutoff::serialize() {
        return {
            {"cutoff", _cutoff},
            {"cutoff_transition_width", _cutoffTransitionWidth}
        };
    }

    double CosCutoff::differentiate(const double r) {
        if (r <= _cutoff - _cutoffTransitionWidth || r >= _cutoff) {
            return 0;
        }
        return -0.5 * M_PI * _cutoffTransitionWidthInverse
                * sin(M_PI*(r - _cutoff + _cutoffTransitionWidth) * _cutoffTransitionWidthInverse) ;
    }

    double CosCutoff::evaluate(const double r) {
        if (r <= _cutoff - _cutoffTransitionWidth) {
            return 1;
        }
        if (r >= _cutoff) {
            return 0;
        }
        return 0.5 * (cos(M_PI * (r - _cutoff + _cutoffTransitionWidth) * _cutoffTransitionWidthInverse) + 1);
    }
}
