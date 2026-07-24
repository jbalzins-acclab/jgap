#ifndef JGAP_DESCRIPTOR_HPP
#define JGAP_DESCRIPTOR_HPP

#include <array>

#include "jgap/core/Real.hpp"

namespace jgap {

    /// Generalized real degrees of freedom of a local atomic descriptor,
    /// without their derivatives.
    ///
    /// \tparam Dim dimensions of the descriptor, not necessarily independent.
    /// \note Just a wrapper around a fixed-sized array, aimed at emphasizing that the array contains descriptor info.
    template<size_t Dim>
    requires (Dim > 0)
    using Descriptor = std::array<Real, Dim>;

}

#endif
