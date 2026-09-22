#ifndef JGAP_TWOBODYDESCRIPTOR_HPP
#define JGAP_TWOBODYDESCRIPTOR_HPP

#include <array>
#include "Descriptor.hpp"
#include "jgap/core/Vector3.hpp"

namespace jgap {

    /// @brief Dim-dimensional 2-body descriptor, and its gradients.
    /// @note Transaltional invariance is assumed, hence grad_r0 = -grad_r1 is not stored explicitly.
    template<size_t Dim>
        requires(Dim > 0)
    struct TwoBodyDescriptor {
        Descriptor<Dim> value{};
        std::array<Vector3, Dim> grad_r1{};

        operator Descriptor<Dim>() { return value; }
        operator Descriptor<Dim>&() { return value; }
        operator const Descriptor<Dim>&() { return value; }
    };

}

#endif
