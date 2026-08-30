#ifndef JGAP_THREEBODYDESCRIPTOR_HPP
#define JGAP_THREEBODYDESCRIPTOR_HPP

#include <array>
#include "Descriptor.hpp"
#include "jgap/core/Vector3.hpp"

namespace jgap {

    /// @brief Dim-dimensional 3-body descriptor and necessary gradients.
    /// @note Assumes translational invariance of the value, hence grad_r0 = -(grad_r1 + grad_r2) is not stored
    /// explicitly.
    template<size_t Dim>
        requires(Dim > 0)
    struct ThreeBodyDescriptor {
        Descriptor<Dim> value{};
        std::array<Vector3, Dim> grad_r1{};
        std::array<Vector3, Dim> grad_r2{};

        operator Descriptor<Dim>() { return value; }
        operator Descriptor<Dim>&() { return value; }
        operator const Descriptor<Dim>&() { return value; }
    };

}

#endif
