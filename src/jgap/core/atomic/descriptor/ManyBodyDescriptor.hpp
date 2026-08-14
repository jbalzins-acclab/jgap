#ifndef JGAP_MANYBODYDESCRIPTOR_HPP
#define JGAP_MANYBODYDESCRIPTOR_HPP

#include <vector>
#include "ThreeBodyDescriptor.hpp"
#include "TwoBodyDescriptor.hpp"
#include "jgap/core/atomic/energy/Virials.hpp"
#include "jgap/core/atomic/geometry/Cluster.hpp"

namespace jgap {

    /// @brief $\vec{Q} = "\sum" \vec{q}_i(n-body-cluster)$, and all possible positional derivatives.
    ///
    /// A local atomic descriptor that may contain information from multiple clusters,
    /// i.e. a sum (or some kind of generalized multiplication) of @ref Descriptor<Dim>
    /// (e.g. a sum of EAM pair-function contributions).
    ///
    /// Stores the descriptor itself, as well as its accumulated derivatives:
    /// forces (negative gradient wrt each atomic position) and @ref Virials.
    ///
    /// @note The number of clusters from which many-body descriptors are constructed
    /// can vary per-atom because of the varying number of neighbours within the cutoff,
    /// so compile-time knowledge of number of per-descriptor derivatives is impossible.
    /// Moreover, a lot of training data usually consists of a very small(2-4) number of atoms,
    /// which, however, corresponds to a rather large number of clusters because of
    /// periodic boundary conditions.
    /// In medium-sized(~100) structures, it is quite usual for around a half of atoms to be
    /// within a cutoff of a single atom.
    /// Hence, accumulating forces/virials is, on average, more efficient than storing a list of
    /// derivatives attributed to each encountered cluster.
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

        /// @brief $\vec{Q} += \vec{q}(r_{ij})$, and update derivatives.
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

        /// @brief $\vec{Q} += \vec{q}(r_{ij}, r_{ik}, r_{jk})$, and update derivatives.
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
