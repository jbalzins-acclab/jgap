#include "core/cutoff/CosCutoff.hpp"

#include <cmath>

#include "utils/Utils.hpp"

namespace jgap {

    CosCutoff::CosCutoff(const double cutoff, const double cutoff_transition_width) {
        cutoff_ = cutoff;
        cutoff_transition_width_ = cutoff_transition_width;
        cutoff_transition_width_inverse_ = 1.0 / cutoff_transition_width_;
    }

    CosCutoff::CosCutoff(const DataNode &params) {
        const auto &cutoffNode = REQUIRE(params, "cutoff");
        cutoff_ = cutoffNode.asDouble();
        if (params.type == DataNode::Type::OBJECT) {
            if (params.contains("cutoff_transition_width")) {
                cutoff_transition_width_ = params.value("cutoff_transition_width", 0.0);
            } else if (params.contains("r_min")) {
                const double rmin = params.value("r_min", 0.0);
                cutoff_transition_width_ = cutoff_ - rmin;
            } else {
                cutoff_transition_width_ = 0.0;
            }
        } else {
            cutoff_transition_width_ = 0.0;
        }
        cutoff_transition_width_inverse_ = cutoff_transition_width_ != 0.0 ? (1.0 / cutoff_transition_width_) : 0.0;
    }

    DataNode CosCutoff::serialize() {
        DataNode obj = DataNode::object();
        auto &m = std::get<std::map<std::string, DataNode>>(obj.value);
        m["cutoff"] = DataNode(cutoff_);
        m["cutoff_transition_width"] = DataNode(cutoff_transition_width_);
        return obj;
    }

    double CosCutoff::differentiate(const double r) {
        if (r <= cutoff_ - cutoff_transition_width_ || r >= cutoff_) {
            return 0;
        }
        return -0.5 * M_PI * cutoff_transition_width_inverse_
                * sin(M_PI*(r - cutoff_ + cutoff_transition_width_) * cutoff_transition_width_inverse_) ;
    }

    double CosCutoff::evaluate(const double r) {
        if (r <= cutoff_ - cutoff_transition_width_) {
            return 1;
        }
        if (r >= cutoff_) {
            return 0;
        }
        return 0.5 * (cos(M_PI * (r - cutoff_ + cutoff_transition_width_) * cutoff_transition_width_inverse_) + 1);
    }
}
