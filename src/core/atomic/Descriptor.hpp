#ifndef JGAP_DESCRIPTOR_HPP
#define JGAP_DESCRIPTOR_HPP

#include "core/atomic/geometry/Cluster.hpp"

#include <array>

namespace jgap {

    template<size_t Dim>
    requires (Dim > 0)
    struct Descriptor {
        std::array<Real, Dim> value{};

        // Modified operator+ to have its own implementation
        Descriptor operator+(const Descriptor& other) const {
            Descriptor result;
            for (size_t i = 0; i < Dim; i++) {
                result.value[i] = value[i] + other.value[i];
            }
            return result;
        }

        // Modified operator+= to use operator+
        Descriptor& operator+=(const Descriptor& other) {
            *this = *this + other;
            return *this;
        }
    };

    template<size_t Dim, size_t DependencyAtoms>
    using DescriptorDerivatives = std::array<std::array<Real, Dim>, Cluster<DependencyAtoms>::NSeparations>;

    template<size_t Dim, size_t DependencyAtoms>
    struct DescriptorAndDerivatives {
        std::array<Real, Dim> value{};
        DescriptorDerivatives<Dim, DependencyAtoms> derivatives{};
    };
}

#endif