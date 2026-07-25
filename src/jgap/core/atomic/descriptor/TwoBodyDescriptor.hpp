#ifndef JGAP_TWOBODYDESCRIPTOR_HPP
#define JGAP_TWOBODYDESCRIPTOR_HPP

#include <array>
#include "Descriptor.hpp"
#include "jgap/core/atomic/geometry/Cluster2.hpp"

namespace jgap {

    /// Derivatives of \ref TwoBodyDescriptor with respect to each separation
    /// in the \ref Cluster2 from which the descriptor was generated.
    template<size_t Dim>
    using TwoBodyDerivatives = std::array<Real, Dim>;

    /// Generalized real degrees of freedom of a local atomic descriptor
    /// derived injectively from a single \ref Cluster2,
    /// as well as their derivatives wrt each separation in the cluster.
    template<size_t Dim>
        requires(Dim > 0)
    struct TwoBodyDescriptor {
        Descriptor<Dim> value{};
        TwoBodyDerivatives<Dim> derivatives{};

        operator Descriptor<Dim>() { return value; }
        operator Descriptor<Dim>&() { return value; }
        operator const Descriptor<Dim>&() { return value; }
    };

} // namespace jgap

#endif
