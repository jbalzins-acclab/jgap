#ifndef JGAP_LATTICE_HPP
#define JGAP_LATTICE_HPP

#include <array>

#include "Vector3.hpp"

namespace jgap {
    struct Lattice {
        Vector3 a, b, c;

        Real volume() const {
            return abs(a.cross(b).dot(c));
        }
    };
}

#endif