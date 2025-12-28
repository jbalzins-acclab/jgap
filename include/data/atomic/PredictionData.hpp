#ifndef JGAP_POTENTIALPREDICTION_HPP
#define JGAP_POTENTIALPREDICTION_HPP

#include "../Vector3.hpp"

namespace jgap {

    struct Predictions {
        std::optional<double> energy;
        std::optional<std::vector<Vector3>> forces;
        std::optional<std::array<Vector3, 3>> virials;

        Predictions operator+(const Predictions& other) const;
    };

    struct Covariance {
        double total;
        std::vector<Vector3> forces;
        std::array<Vector3, 3> virials;
    };
}

#endif