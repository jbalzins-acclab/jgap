#ifndef JGAP_MANYBODYDESCRIPTOR_HPP
#define JGAP_MANYBODYDESCRIPTOR_HPP

#include <vector>
#include "ThreeBodyDescriptor.hpp"
#include "TwoBodyDescriptor.hpp"
#include "jgap/core/atomic/energy/Virials.hpp"
#include "jgap/core/atomic/geometry/Cluster.hpp"

namespace jgap {

    /// Generalized real degrees of freedom of a local atomic descriptor
    /// that contains information from a multitude of \ref Clusters,
    /// i.e. ~ a sum/some kind of multiplication of \ref NBodyDescriptors
    /// (e.g. a sum of EAM pair-function contributions).
    ///
    /// Stores the descriptor itself, as well as its accumulated forces per atom and virials.
    ///
    /// \tparam Dim dimensions of the descriptor
    ///
    /// \note The number of clusters from which many-body descriptors are constructed
    /// can vary per-atom because of the varying number of neighbours within the cutoff,
    /// so compile-time knowledge of per-cluster derivative number is impossible.
    /// Moreover, a lot of training data usually consists of a very small(2-4) number of atoms,
    /// which, however, corresponds to a rather large number of clusters because of
    /// the periodic boundary conditions.
    /// In medium-sized(~100) structures, it is quite usual for around a half of atoms to be
    /// within a cutoff of a single atom.
    /// Hence, accumulating forces/virials is in general more effective than storing a list of
    /// DescriptorDerivatives per cluster.
    template<size_t Dim>
        requires(Dim > 0)
    struct ManyBodyDescriptor {
        Descriptor<Dim> value{};
        std::array<Virials, Dim> virials{};
        std::vector<std::array<Vector3, Dim>> forces{};

        ManyBodyDescriptor(size_t n_atoms) : forces(n_atoms) {}

        operator Descriptor<Dim>() { return value; }
        operator Descriptor<Dim>&() { return value; }
        operator const Descriptor<Dim>&() { return value; }

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

        ManyBodyDescriptor operator+(const ManyBodyDescriptor& other) const {
            ManyBodyDescriptor result = *this;
            result += other;
            return result;
        }

        ManyBodyDescriptor& operator*=(Real scalar) {
            for (size_t d = 0; d < Dim; d++) {
                value[d] *= scalar;
                virials[d] *= scalar;
            }
            for (size_t i = 0; i < forces.size(); i++) {
                for (size_t d = 0; d < Dim; d++) {
                    forces[i][d] *= scalar;
                }
            }
            return *this;
        }

        ManyBodyDescriptor operator*(Real scalar) const {
            ManyBodyDescriptor result = *this;
            result *= scalar;
            return result;
        }

        /// \usage
        void add(const Cluster2& cluster, const TwoBodyDescriptor<Dim>& contribution) {
            for (size_t d = 0; d < Dim; d++) {
                value[d] += contribution.value[d];
            }

            const auto& separation = cluster.separation01;
            const auto& derivs = contribution.derivatives;

            for (size_t dim = 0; dim < Dim; dim++) {
                Vector3 force_contrib = separation.direction * derivs[dim];
                forces[cluster.idx0][dim] += force_contrib;
                forces[cluster.idx1][dim] -= force_contrib;

                virials[dim] += separation.virials() * derivs[dim];
            }
        }

        void add(const Cluster3& cluster, const ThreeBodyDescriptor<Dim>& contribution) {
            for (size_t d = 0; d < Dim; d++) {
                value[d] += contribution.value[d];
            }

            for (size_t i = 0; i < 3; i++) {
                for (size_t j = i + 1; j < 3; j++) {
                    const auto sep_idx = flattenedIndex(i, j);
                    const auto& separation = cluster.separation(sep_idx);
                    const auto& derivs = contribution.derivatives[sep_idx];

                    for (size_t dim = 0; dim < Dim; dim++) {
                        Vector3 force_contrib = separation.direction * derivs[dim];
                        forces[cluster.atom_indexes[i]][dim] += force_contrib;
                        forces[cluster.atom_indexes[j]][dim] -= force_contrib;

                        virials[dim] += separation.virials() * derivs[dim];
                    }
                }
            }
        }
    };

}

#endif
