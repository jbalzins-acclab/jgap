#ifndef JGAP_LATTICE_HPP
#define JGAP_LATTICE_HPP

#include <array>

#include "../../Vector3.hpp"

namespace jgap {
    /// @brief Periodic lattice vectors.
    /// @note A vector is treated as a dummy variable if PBC is set to false in the respective dimension in \ref Atoms.
    struct Lattice {
        Vector3 a, b, c;

        Real volume() const {
            return abs(a.cross(b).dot(c));
        }
    };
}

#endif