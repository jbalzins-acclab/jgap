#ifndef JGAP_SQUAREDEXPKERNEL_HPP
#define JGAP_SQUAREDEXPKERNEL_HPP

#include "Kernel.hpp"
#include <cmath>

namespace jgap {

    template<size_t MainDimensions, size_t CutoffDimensions, size_t DescriptorDeps>
    class SquaredExpKernel : public Kernel<MainDimensions+CutoffDimensions, DescriptorDeps> {
    public:
        static constexpr size_t TotalDimensions = MainDimensions + CutoffDimensions;

        using KernelValueAndGradient = Kernel<MainDimensions+CutoffDimensions, DescriptorDeps>::KernelValueAndGradient;

        SquaredExpKernel() = default;

        SquaredExpKernel(const Real energy_scale,
                         const std::array<Real, MainDimensions>& length_scales)
        {
            prefactor = energy_scale * energy_scale;
            for (size_t dim = 0; dim < MainDimensions; dim++) {
                inverse_length_scales_squared[dim] = 1.0 / (length_scales[dim] * length_scales[dim]);
            }
        }

        Real value(const Kernel<TotalDimensions, DescriptorDeps>::TValue &q1,
                   const Kernel<TotalDimensions, DescriptorDeps>::TValue &q2) const override {

            Real exp_argument = 0;
            for (size_t dim = 0; dim < MainDimensions; dim++) {
                Real diff = q1[dim] - q2[dim];
                exp_argument += diff * diff * inverse_length_scales_squared[dim];
            }
            Real val = prefactor * std::exp(-0.5 * exp_argument);

            for (size_t dim = MainDimensions; dim < TotalDimensions; dim++) {
                val *= q1[dim] * q2[dim];
            }

            return val;
        }

        KernelValueAndGradient valueAndGradient(
            const Kernel<TotalDimensions, DescriptorDeps>::TValue &sparse_point,
            const Kernel<TotalDimensions, DescriptorDeps>::TValue &q) const override {

            Real exp_argument = 0;
            for (size_t dim = 0; dim < MainDimensions; dim++) {
                Real diff = q[dim] - sparse_point[dim];
                exp_argument += diff * diff * inverse_length_scales_squared[dim];
            }
            Real val = prefactor * std::exp(-0.5 * exp_argument);

            for (size_t dim = MainDimensions; dim < TotalDimensions; dim++) {
                val *= sparse_point[dim];
            }

            std::array<Real, TotalDimensions> gradient{};
            for (size_t dim = MainDimensions; dim < TotalDimensions; dim++) {
                gradient[dim] = val;
            }

            for (size_t dim = MainDimensions; dim < TotalDimensions; dim++) {
                val *= q[dim];

                for (size_t j = MainDimensions; j < dim; j++) {
                    gradient[j] *= q[dim];
                }

                for (size_t j = dim + 1; j < TotalDimensions; j++) {
                    gradient[j] *= q[dim];
                }
            }

            for (size_t dim = 0; dim < MainDimensions; dim++) {
                gradient[dim] = val * (sparse_point[dim] - q[dim]) * inverse_length_scales_squared[dim];
            }

            return {
                .value = val,
                .gradient = gradient
            };
        }

    private:
        Real prefactor{};
        std::array<Real, MainDimensions> inverse_length_scales_squared{};
    };
}

#endif