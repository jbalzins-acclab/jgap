#ifndef JGAP_THREEBODYDESCRIPTOR_HPP
#define JGAP_THREEBODYDESCRIPTOR_HPP

#include "Descriptor.hpp"
#include "jgap/core/atomic/geometry/Cluster3.hpp"
#include <array>

namespace jgap {

    /// @brief Derivatives of Descriptor<Dim>(r_ij, r_ik, r_jk) wrt (r_ij, r_ik, r_jk).
    template<size_t Dim>
    using ThreeBodyDerivatives = std::array<std::array<Real, Dim>, 3>;

    /// @brief Dim-dimensional 3-body descriptor, and derivatives wrt (r_ij, r_ik, r_jk) in each dimension.
    template<size_t Dim>
        requires(Dim > 0)
    struct ThreeBodyDescriptor {
        Descriptor<Dim> value{};
        ThreeBodyDerivatives<Dim> derivatives{};

        operator Descriptor<Dim>() { return value; }
        operator Descriptor<Dim> &() { return value; }
        operator const Descriptor<Dim> &() { return value; }
    };

} // namespace jgap

#endif
