#ifndef JGAP_POTENTIALPREDICTION_HPP
#define JGAP_POTENTIALPREDICTION_HPP

#include "Vector3.hpp"

namespace jgap {

    struct Predictions {
        optional<double> energy;
        optional<vector<Vector3>> forces;
        optional<array<Vector3, 3>> virials;

        Predictions operator+(const Predictions& other) const;
    };

    struct Covariance {
        double total;
        vector<Vector3> forces;
        array<Vector3, 3> virials;
    };
}

#endif