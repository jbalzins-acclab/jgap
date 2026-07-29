#ifndef JGAP_MANYBODYDESCRIPTORS_HPP
#define JGAP_MANYBODYDESCRIPTORS_HPP

#include <cassert>
#include <vector>
#include "ManyBodyDescriptor.hpp"
#include "ThreeBodyDescriptor.hpp"
#include "TwoBodyDescriptor.hpp"
#include "jgap/core/Vector3.hpp"
#include "jgap/core/atomic/geometry/Cluster.hpp"

namespace jgap {

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

        std::array<Vector3, Dim>& forces_at(size_t desc_idx, size_t atom_idx) {
            return forces[desc_idx * n_atoms + atom_idx];
        }

        const std::array<Vector3, Dim>& forces_at(size_t desc_idx, size_t atom_idx) const {
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
                for (size_t d = 0; d < Dim; d++) {
                    values[i][d] *= scalar;
                    virials[i][d] *= scalar;
                }
            }
            for (size_t i = 0; i < forces.size(); i++) {
                for (size_t d = 0; d < Dim; d++) {
                    forces[i][d] *= scalar;
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
                for (size_t d = 0; d < Dim; d++) {
                    total.value[d] += values[i][d] * coeffs[i];
                    total.virials[d] += virials[i][d] * coeffs[i];
                }

                for (size_t j = 0; j < n_atoms; j++) {
                    for (size_t d = 0; d < Dim; d++) {
                        total.forces[j][d] += forces[i * n_atoms + j][d] * coeffs[i];
                    }
                }
            }
            return total;
        }

        void add(size_t descriptor_index, const Cluster2& cluster, const TwoBodyDescriptor<Dim>& contribution) {
            for (size_t d = 0; d < Dim; d++) {
                values[descriptor_index][d] += contribution.value[d];
            }

            if (contribution.derivatives.empty()) {
                return;
            }

            const auto& separation = cluster.separation01;
            const auto& derivs = contribution.derivatives;

            for (size_t dim = 0; dim < Dim; dim++) {
                Vector3 force_contrib = separation.direction * derivs[dim];
                forces[descriptor_index * n_atoms + cluster.idx0][dim] += force_contrib;
                forces[descriptor_index * n_atoms + cluster.idx1][dim] -= force_contrib;

                virials[descriptor_index][dim] += separation.virials() * derivs[dim];
            }
        }

        void add(size_t descriptor_index, const Cluster3& cluster, const ThreeBodyDescriptor<Dim>& contribution) {
            for (size_t d = 0; d < Dim; d++) {
                values[descriptor_index][d] += contribution.value[d];
            }

            if (contribution.derivatives.empty()) {
                return;
            }

            for (size_t i = 0; i < 3; i++) {
                for (size_t j = i + 1; j < 3; j++) {
                    const auto sep_idx = flattenedIndex(i, j);
                    const auto& separation = cluster.separation(sep_idx);
                    const auto& derivs = contribution.derivatives[sep_idx];

                    for (size_t dim = 0; dim < Dim; dim++) {
                        Vector3 force_contrib = separation.direction * derivs[dim];
                        forces[descriptor_index * n_atoms + cluster.atom_indexes[i]][dim] += force_contrib;
                        forces[descriptor_index * n_atoms + cluster.atom_indexes[j]][dim] -= force_contrib;

                        virials[descriptor_index][dim] += separation.virials() * derivs[dim];
                    }
                }
            }
        }
    };

}

#endif
