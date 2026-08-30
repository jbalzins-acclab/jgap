#ifndef JGAP_MANYBODYDESCRIPTORS_HPP
#define JGAP_MANYBODYDESCRIPTORS_HPP

#include <cassert>
#include <vector>
#include "ManyBodyDescriptor.hpp"
#include "ThreeBodyDescriptor.hpp"
#include "TwoBodyDescriptor.hpp"
#include "jgap/core/Vector3.hpp"
#include "jgap/core/atomic/geometry/Cluster2.hpp"
#include "jgap/core/atomic/geometry/Cluster3.hpp"
#include "jgap/utils/Utils.hpp"

namespace jgap {

    /// @brief An array of @ref ManyBodyDescriptor with separate descriptor properties stored in a continuous block.
    /// @note Helps to avoid strong a vector of pointers (separate force vectors).
    template<size_t Dim>
        requires(Dim > 0)
    struct ManyBodyDescriptors {
        size_t n_atoms = 0;
        std::vector<Descriptor<Dim>> values{};
        std::vector<std::array<Virials, Dim>> virials{};
        std::vector<std::array<Vector3, Dim>> forces{};

        ManyBodyDescriptors() = default;

        ManyBodyDescriptors(size_t n_descriptors, size_t n_atoms) :
            n_atoms(n_atoms), values(n_descriptors), virials(n_descriptors), forces(n_descriptors * n_atoms) {}

        std::array<Vector3, Dim>& force(size_t desc_idx, size_t atom_idx) {
            return forces[desc_idx * n_atoms + atom_idx];
        }

        const std::array<Vector3, Dim>& force(size_t desc_idx, size_t atom_idx) const {
            return forces[desc_idx * n_atoms + atom_idx];
        }

        ManyBodyDescriptors& operator+=(const ManyBodyDescriptors& other) {
            assert(n_atoms == other.n_atoms);
            assert(values.size() == other.values.size());

            for (size_t i = 0; i < values.size(); i++) {
                for (size_t d = 0; d < Dim; d++) {
                    values[i][d] += other.values[i][d];
                    virials[i][d] += other.virials[i][d];
                }
            }
            for (size_t i = 0; i < forces.size(); i++) {
                for (size_t d = 0; d < Dim; d++) {
                    forces[i][d] += other.forces[i][d];
                }
            }
            return *this;
        }

        ManyBodyDescriptors operator+(const ManyBodyDescriptors& other) const {
            ManyBodyDescriptors result = *this;
            result += other;
            return result;
        }

        ManyBodyDescriptors& operator*=(Real scalar) {
            for (size_t i = 0; i < values.size(); i++) {
                for (size_t dim = 0; dim < Dim; dim++) {
                    values[i][dim] *= scalar;
                    virials[i][dim] *= scalar;
                }
            }
            for (size_t i = 0; i < forces.size(); i++) {
                for (size_t dim = 0; dim < Dim; dim++) {
                    forces[i][dim] *= scalar;
                }
            }
            return *this;
        }

        ManyBodyDescriptors operator*(Real scalar) const {
            ManyBodyDescriptors result = *this;
            result *= scalar;
            return result;
        }

        ManyBodyDescriptor<Dim> sum(const std::vector<Real>& coeffs) const {
            assert(coeffs.size() == values.size());
            ManyBodyDescriptor<Dim> total(n_atoms);

            for (size_t i = 0; i < values.size(); i++) {
                for (size_t dim = 0; dim < Dim; dim++) {
                    total.value[dim] += values[i][dim] * coeffs[i];
                    total.virials[dim] += virials[i][dim] * coeffs[i];
                }

                for (size_t j = 0; j < n_atoms; j++) {
                    for (size_t dim = 0; dim < Dim; dim++) {
                        total.forces[j][dim] += forces[i * n_atoms + j][dim] * coeffs[i];
                    }
                }
            }
            return total;
        }

        /// @brief $\vec{Q}_{descriptor-index} += \vec{q}(\vec{r}_{ij})$, and update derivatives.
        void add(size_t descriptor_index, const Cluster2& cluster, const TwoBodyDescriptor<Dim>& contribution) {
            const auto r01 = cluster.separation01.vec();

            for (size_t dim = 0; dim < Dim; dim++) {
                values[descriptor_index][dim] += contribution.value[dim];

                forces[descriptor_index * n_atoms + cluster.idx1][dim] -= contribution.grad_r1[dim];
                forces[descriptor_index * n_atoms + cluster.idx0][dim] += contribution.grad_r1[dim];

                virials[descriptor_index][dim] += Virials::dyadic(r01, -contribution.grad_r1[dim]);
            }
        }

        /// @brief $\vec{Q}_{descriptor-index} += \vec{q}(\vec{r}_{ij}, \vec{r}_{ik}, \vec{r}_{jk})$, and update
        /// derivatives.
        void add(size_t descriptor_index, const Cluster3& cluster, const ThreeBodyDescriptor<Dim>& contribution) {
            const auto r01 = cluster.separation01.vec();
            const auto r02 = cluster.separation02.vec();

            for (size_t dim = 0; dim < Dim; dim++) {
                values[descriptor_index][dim] += contribution.value[dim];

                forces[descriptor_index * n_atoms + cluster.idx1][dim] -= contribution.grad_r1[dim];
                forces[descriptor_index * n_atoms + cluster.idx2][dim] -= contribution.grad_r2[dim];
                forces[descriptor_index * n_atoms + cluster.idx0][dim] +=
                    (contribution.grad_r1[dim] + contribution.grad_r2[dim]);

                virials[descriptor_index][dim] += Virials::dyadic(r01, -contribution.grad_r1[dim]);
                virials[descriptor_index][dim] += Virials::dyadic(r02, -contribution.grad_r2[dim]);
            }
        }
    };

}

#endif
