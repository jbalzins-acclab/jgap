#ifndef JGAP_REAL_HPP
#define JGAP_REAL_HPP

#include <type_traits>

namespace jgap {
    #ifndef FPTYPE
        using Real = double;
    #else
        using Real = FPTYPE;
    #endif
}

// Just in case
static_assert(std::is_floating_point_v<jgap::Real>, "Real must be a floating point type.");

#endif