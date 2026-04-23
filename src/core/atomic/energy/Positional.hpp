#ifndef JGAP_PREDICTIONS_HPP
#define JGAP_PREDICTIONS_HPP

#include "../geometry/Vector3.hpp"
#include <vector>

#include "Virials.hpp"

namespace jgap {
    struct Positional {
        // itself
        Real value;
        // derivatives
        Virials virials;
        std::vector<Vector3> forces;

        Positional operator+(const Positional& other) const;
        Positional& operator+=(const Positional& other);
        Positional operator*(double scalar) const;
        Positional& operator*=(double scalar);

        Positional operator*(const Positional& other) const;
    };
}

#endif
