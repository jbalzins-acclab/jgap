#ifndef JGAP_REAL_HPP
#define JGAP_REAL_HPP

#include <type_traits>

namespace jgap {
    /// @brief Type Alias for the floating point number type.
    /// @note Meant to SIMPLIFY transition from double to float if such is needed in the future,
    /// it is not guaranteed to work by simple switching.
    /// @note Initial testing of floats showed ~20% speedup on a small single element potential,
    /// returning NaNs as a result. It's probably not worth any further testing,
    /// except using floats to store the covariance matrix but using doubles for operations when finding least squares
    /// (like in QRGapFit) may work.
    using Real = double;

    inline namespace literals {
        constexpr Real operator""_r(long double val) { return static_cast<Real>(val); }
        constexpr Real operator""_r(unsigned long long val) { return static_cast<Real>(val); }
    }
}

// Just in case
static_assert(std::is_floating_point_v<jgap::Real>, "Real must be a floating point type.");

#endif
