#include "core/cutoff/CosCutoff.hpp"

#include <cmath>

#include "utils/Utils.hpp"

namespace jgap {

    CosCutoff::CosCutoff(const Real cutoff, const Real cutoff_transition_width)
        : cutoff(cutoff), cutoff_transition_width(cutoff_transition_width)
    {
        cutoff_transition_width_inverse = 1.0 / cutoff_transition_width;
    }

    /*CosCutoff::CosCutoff(const DataNode &params) {
        const auto &cutoffNode = REQUIRE(params, "cutoff");
        cutoff = cutoffNode.asReal();
        if (params.type == DataNode::Type::OBJECT) {
            if (params.contains("cutoff_transition_width")) {
                cutoff_transition_width = params.value("cutoff_transition_width", 0.0);
            } else if (params.contains("r_min")) {
                const Real rmin = params.value("r_min", 0.0);
                cutoff_transition_width = cutoff - rmin;
            } else {
                cutoff_transition_width = 0.0;
            }
        } else {
            cutoff_transition_width = 0.0;
        }
        cutoff_transition_width_inverse = cutoff_transition_width != 0.0 ? (1.0 / cutoff_transition_width) : 0.0;
    }

    DataNode CosCutoff::serialize() {
        DataNode obj = DataNode::object();
        auto &m = std::get<std::map<std::string, DataNode>>(obj.value);
        m["cutoff"] = DataNode(cutoff);
        m["cutoff_transition_width"] = DataNode(cutoff_transition_width);
        return obj;
    }*/

    Real CosCutoff::differentiate(const Real r) const {
        if (r <= cutoff - cutoff_transition_width || r >= cutoff) {
            return 0;
        }
        return - static_cast<Real>(0.5) * M_PI * cutoff_transition_width_inverse
                * sin(M_PI*(r - cutoff + cutoff_transition_width) * cutoff_transition_width_inverse) ;
    }

    Real CosCutoff::evaluate(const Real r) const {
        if (r <= cutoff - cutoff_transition_width) {
            return 1;
        }
        if (r >= cutoff) {
            return 0;
        }
        return static_cast<Real>(0.5) * (
            cos(M_PI * (r - cutoff + cutoff_transition_width) * cutoff_transition_width_inverse)
            + static_cast<Real>(1.0)
            );
    }
}
