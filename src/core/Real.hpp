#ifndef JGAP_REAL_HPP
#define JGAP_REAL_HPP

#include <type_traits>

namespace jgap {
    /**
     * Type Alias for the floating point number type.
     * Meant to SIMPLIFY transition from double to float if such is needed in the future.
     * However, further adjustments in the code will be needed,
     * as I found putting static_cast everywhere quite annoying,
     * yet the testing I ran showed not much benefit from such transition with a standard QR fit
     * (no notable speedup, yet nonsense coefficients).
     */
    using Real = double;
}

// Just in case
static_assert(std::is_floating_point_v<jgap::Real>, "Real must be a floating point type.");

#endif