#ifndef JGAP_THREEBODYDESCRIPTOR_HPP
#define JGAP_THREEBODYDESCRIPTOR_HPP

#include "Descriptor.hpp"
#include "jgap/core/atomic/geometry/Cluster3.hpp"
#include <array>

namespace jgap {

    /// Derivatives of \ref ThreeBodyDescriptor with respect to each separation
    /// in the \ref Cluster3 from which the descriptor was generated.
    template<size_t Dim>
    using ThreeBodyDerivatives = std::array<std::array<Real, Dim>, 3>;

    /// Generalized real degrees of freedom of a local atomic descriptor
    /// derived injectively from a single \ref Cluster3,
    /// as well as their derivatives wrt each separation in the cluster.
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
