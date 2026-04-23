#ifndef JGAP_CUTOFFRANGES_HPP
#define JGAP_CUTOFFRANGES_HPP

#include <optional>

#include "utils/Utils.hpp"

namespace jgap {
    struct CutoffRanges {
        std::optional<double> two_body;
        std::optional<double> three_body;
        std::optional<double> min_eam_density, max_eam_density;

        double maxOverall() const {
            double result = 0.0;
            if (two_body.has_value()) result = two_body.value();
            if (three_body.has_value()) result = three_body.value();
            return result;
        }

        CutoffRanges operator+(const CutoffRanges& other) const {
            auto res = *this;
            if (other.two_body.has_value()) {
                res.two_body = std::max(res.two_body.value_or(0.0), other.two_body.value());
            }
            if (other.three_body.has_value()) {
                res.three_body = std::max(res.three_body.value_or(0.0), other.three_body.value());
            }
            if (other.max_eam_density.has_value()) {
                res.max_eam_density = std::max(res.max_eam_density.value_or(0.0), other.max_eam_density.value());
            }
            if (other.min_eam_density.has_value()) {
                res.min_eam_density = std::min(res.min_eam_density.value_or(0.0), other.min_eam_density.value());
            }
            return res;
        }

        CutoffRanges& operator+=(const CutoffRanges& other) {
            *this = *this + other;
            return *this;
        }
    };
}

#endif