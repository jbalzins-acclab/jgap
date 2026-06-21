#ifndef JGAP_DESCRIPTOR_HPP
#define JGAP_DESCRIPTOR_HPP

#include "core/atomic/geometry/Cluster.hpp"

#include <array>

namespace jgap {

    /**
     * Generalized real degrees of freedom of a local atomic descriptor,
     * without their derivatives.
     *
     * @tparam Dim dimensions of the descriptor, not necessarily independent.
     * @note Just a wrapper around a fixed-sized array,
     * aimed at emphasizing that the array contains descriptor info
     * and to keep the same style with {@link NBodyDescriptor} and {@link ManyBodyDescriptor}.
     */
    template<size_t Dim>
    requires (Dim > 0)
    struct Descriptor {
        std::array<Real, Dim> value{};
    };

    /**
     * Derivatives of {@link NBodyDescriptor} with respect to each separation
     * in the {@link Cluster} from which the descriptor was generated.
     *
     * Stored in format: derivative[wrt n-th separation in cluster][of m-th descriptor's dimension]
     *
     * @tparam Dim dimensions of the descriptor
     * @tparam ClusterSize size of the cluster
     */
    template<size_t Dim, size_t ClusterSize>
    using DescriptorDerivatives = std::array<std::array<Real, Dim>, Cluster<ClusterSize>::NSeparations>;

    /**
     * Generalized real degrees of freedom of a local atomic descriptor
     * derived injectively from a single {@link Cluster},
     * as well as their derivatives wrt each separation in the cluster.
     *
     * @tparam Dim dimensions of the descriptor
     * @tparam ClusterSize size of the cluster
     *
     * @note Layout is intended to simplify conversion to {@link Descriptor}.
     */
    template<size_t Dim, size_t ClusterSize>
    requires (Dim > 0 && ClusterSize > 0)
    struct NBodyDescriptor {
        std::array<Real, Dim> value{};
        DescriptorDerivatives<Dim, ClusterSize> derivatives{};
    };

    /**
     * Generalized real degrees of freedom of a local atomic descriptor
     * that contains information from a multitude of {@link Cluster}s,
     * i.e. ~ a sum/some kind of multiplication of {@link NBodyDescriptor}s
     * like a sum of EAM pair-function contributions.
     *
     * Stores the descriptor itself, as well as its accumulated forces per atom and virials.
     *
     * @tparam Dim dimensions of the descriptor
     *
     * @note The number of clusters from which many-body descriptors are constructors
     * can vary per-atom because of the varying number of neighbours within the cutoff,
     * so compile-time knowledge of per-cluster derivative number is impossible.
     * Moreover, a lot of training data usually consists of a very small(2-4) number of atoms,
     * which, however, corresponds to a rather large number of clusters because of
     * the periodic boundary conditions.
     * In medium-sized(~100) structures, it is quite usual for around a half of atoms to be
     * within a cutoff of a single atom.
     * Hence, accumulating forces/virials is in general more effective than storing a list of
     * {@link DescriptorDerivatives} per cluster.
     */
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