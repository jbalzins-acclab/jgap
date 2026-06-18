#ifndef JGAP_DESCRIPTOR_HPP
#define JGAP_DESCRIPTOR_HPP

#include "core/atomic/geometry/Cluster.hpp"

#include <array>

namespace jgap {

    template<size_t Dim>
    requires (Dim > 0)
    struct Descriptor {
        std::array<Real, Dim> value{};
    };

    template<size_t Dim, size_t ClusterSize>
    using DescriptorDerivatives = std::array<std::array<Real, Dim>, Cluster<ClusterSize>::NSeparations>;

    template<size_t Dim, size_t ClusterSize>
    requires (Dim > 0 && ClusterSize > 0)
    struct NBodyDescriptor {
        std::array<Real, Dim> value{};
        DescriptorDerivatives<Dim, ClusterSize> derivatives{};
    };

    template<size_t Dim>
    requires (Dim > 0)
    struct ManyBodyDescriptor {
        std::array<Real, Dim> value{};
        std::array<Virials, Dim> virials{};
        std::vector<std::array<Vector3, Dim>> forces{};

        ManyBodyDescriptor() = default;
        ManyBodyDescriptor(size_t n_atoms) : forces(n_atoms) {}

        ManyBodyDescriptor& operator+=(const ManyBodyDescriptor& other) {
            for (size_t d = 0; d < Dim; d++) {
                value[d] += other.value[d];
                virials[d] += other.virials[d];
            }
            for (size_t i = 0; i < forces.size(); i++) {
                for (size_t d = 0; d < Dim; d++) {
                    forces[i][d] += other.forces[i][d];
                }
            }
            return *this;
        }
    };
}

#endif