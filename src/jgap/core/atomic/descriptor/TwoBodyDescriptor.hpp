#ifndef JGAP_TWOBODYDESCRIPTOR_HPP
#define JGAP_TWOBODYDESCRIPTOR_HPP

#include <array>
#include "Descriptor.hpp"
#include "jgap/core/atomic/geometry/Cluster2.hpp"

namespace jgap {

    /// @brief Derivatives of Descriptor<Dim>(r_ij) wrt r_ij.
    template<size_t Dim>
    using TwoBodyDerivatives = std::array<Real, Dim>;

    /// @brief Dim-dimensional 2-body descriptor, and derivatives wrt r_ij in each dimension.
    template<size_t Dim>
        requires(Dim > 0)
    struct TwoBodyDescriptor {
        Descriptor<Dim> value{};
        TwoBodyDerivatives<Dim> derivatives{};

        operator Descriptor<Dim>() { return value; }
        operator Descriptor<Dim>&() { return value; }
        operator const Descriptor<Dim>&() { return value; }
    };

}

#endif
