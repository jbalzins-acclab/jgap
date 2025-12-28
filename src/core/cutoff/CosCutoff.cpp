#include "core/cutoff/CosCutoff.hpp"

#include <cmath>

#include "utils/Utils.hpp"

namespace jgap {

    CosCutoff::CosCutoff(const double cutoff, const double cutoffTransitionWidth) {
        _cutoff = cutoff;
        _cutoffTransitionWidth = cutoffTransitionWidth;
        _cutoffTransitionWidthInverse = 1.0 / _cutoffTransitionWidth;
    }

    CosCutoff::CosCutoff(const DataNode &params) {
        const auto &cutoffNode = require(params, "cutoff");
        _cutoff = cutoffNode.asDouble();
        if (params.type == DataNode::Type::OBJECT) {
            if (params.contains("cutoff_transition_width")) {
                _cutoffTransitionWidth = params.value("cutoff_transition_width", 0.0);
            } else if (params.contains("r_min")) {
                const double rmin = params.value("r_min", 0.0);
                _cutoffTransitionWidth = _cutoff - rmin;
            } else {
                _cutoffTransitionWidth = 0.0;
            }
        } else {
            _cutoffTransitionWidth = 0.0;
        }
        _cutoffTransitionWidthInverse = _cutoffTransitionWidth != 0.0 ? (1.0 / _cutoffTransitionWidth) : 0.0;
    }

    DataNode CosCutoff::serialize() {
        DataNode obj = DataNode::object();
        auto &m = std::get<std::map<std::string, DataNode>>(obj.value);
        m["cutoff"] = DataNode(_cutoff);
        m["cutoff_transition_width"] = DataNode(_cutoffTransitionWidth);
        return obj;
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
