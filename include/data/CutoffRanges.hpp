#ifndef JGAP_RANGE_HPP
#define JGAP_RANGE_HPP

#include <optional>

#include "utils/Utils.hpp"


namespace jgap {
    struct CutoffRanges {
        optional<double> twoBody;
        optional<double> threeBody;
        optional<double> minEam, maxEam;

        double maxOverall() const {
            double result = 0;
            if (twoBody.has_value()) result = twoBody.value();
            if (threeBody.has_value()) result = threeBody.value();
            return result;
        }

        CutoffRanges operator+(const CutoffRanges& other) const {
            auto res = *this;
            if (other.twoBody.has_value()) {
                res.twoBody = max(res.twoBody.value_or(0.0), other.twoBody.value());
            }
            if (other.threeBody.has_value()) {
                res.threeBody = max(res.threeBody.value_or(0.0), other.threeBody.value());
            }
            if (other.maxEam.has_value()) {
                res.maxEam = max(res.maxEam.value_or(0.0), other.maxEam.value());
            }
            if (other.minEam.has_value()) {
                res.minEam = min(res.minEam.value_or(0.0), other.minEam.value());
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