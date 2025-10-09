#include "core/cutoff/CosCutoff.hpp"

#include <cmath>

namespace jgap {

    CosCutoff::CosCutoff(nlohmann::json params) {
        _cutoff = params["cutoff"];
        if (params.contains("cutoff_transition_width")) {
            _cutoffTransitionWidth = params["cutoff_transition_width"];
        } else {
            _cutoffTransitionWidth = _cutoff - params["r_min"].get<double>();
        }
        _cutoffTransitionWidthInverse = (1.0 / _cutoffTransitionWidth);
    }

    CosCutoff::CosCutoff(const double cutoff, const double rMin) {
        _cutoff = cutoff;
        _cutoffTransitionWidth = cutoff - rMin;
        _cutoffTransitionWidthInverse = (1.0 / _cutoffTransitionWidth);
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
