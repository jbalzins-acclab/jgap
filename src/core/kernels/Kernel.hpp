#ifndef JGAP_KERNEL_HPP
#define JGAP_KERNEL_HPP

#include "../atomic/geometry/Vector3.hpp"
#include <optional>
#include <string>
#include <memory>
#include <array>
#include <vector>

#include "../atomic/Descriptor.hpp"
#include "core/atomic/energy/AtomicQuantity.hpp"

namespace jgap {

    template<size_t DescriptorDim, size_t DescriptorDependencies>
    class Kernel {
    public:
        constexpr static size_t Dim = DescriptorDim;
        constexpr static size_t Dependencies = DescriptorDependencies;

        using TDescriptor = Descriptor<DescriptorDim, DescriptorDependencies, CalculationType::ValueOnly>;
        using TDescWithGradients = Descriptor<DescriptorDim, DescriptorDependencies, CalculationType::WithGradients>;
        using TValue = typeof(TDescriptor::value);

        struct KernelValueAndGradient {
            Real value;
            std::array<Real, Dim> gradient;
        };

        virtual ~Kernel() = default;

        virtual Real value(const TValue& q1, const TValue& q2) const = 0;
        virtual KernelValueAndGradient valueAndGradient(const TValue& sparse_point, const TValue& q) const = 0;

        Real covariance(const TValue& sparse_point, const std::vector<TDescriptor>& descriptors) const;
        void covariance(AtomicQuantity& added_to, const TValue& sparse_point,
                        const std::vector<TDescWithGradients>& descriptors) const;
    };

    template<size_t DescriptorDim, size_t DescriptorDependencies>
    void Kernel<DescriptorDim, DescriptorDependencies>::covariance(AtomicQuantity &added_to,
        const TValue &sparse_point, const std::vector<TDescWithGradients> &descriptors) const {

        Real& energy = added_to.value;
        Virials& virials = added_to.virials;
        std::vector<Vector3>& forces = added_to.forces;

        for (const auto &descriptor: descriptors) {

            // g[i] = dU/dq_i
            const auto [val, gradK_wrt_q] = valueAndGradient(sparse_point, descriptor.value);

            energy += val;

            for (size_t dim = 0; dim < DescriptorDim; dim++) {
                virials += descriptor.virials[dim] * gradK_wrt_q[dim];
            }

            if constexpr (DescriptorDependencies == GradientData::UNKNOWN_DEPENDENCIES) {
                // Descriptor contains: std::array<std::vector<Vector3>, Dim> gradients;
                for (int i = 0; i < forces.size(); i++) {
                    for (size_t dim = 0; dim < DescriptorDim; dim++) {
                        forces[i] -= descriptor.gradients[dim][i] * gradK_wrt_q[dim];
                    }
                }
            } else {
                // Descriptor contains: std::array<std::array<GradientData, DependencyAtoms>, Dim> gradients{};
                for (size_t dim = 0; dim < DescriptorDim; dim++) {
                    for (size_t i = 0; i < DescriptorDependencies; i++) {
                        const GradientData& gradQdim_wrt_ri = descriptor.gradients[dim][i];
                        forces[gradQdim_wrt_ri.wrt_atom_index] -= gradQdim_wrt_ri.value * gradK_wrt_q[dim];
                    }
                }
            }
        }
    }


    template<size_t DescriptorDim, size_t DescriptorDependencies>
    Real Kernel<DescriptorDim, DescriptorDependencies>::covariance(
        const TValue &sparse_point, const std::vector<TDescriptor> &descriptors) const {

        Real energy = 0;

        for (const auto &descriptor: descriptors) {
            energy += value(sparse_point, descriptor.value);
        }

        return energy;
    }
}

#endif
