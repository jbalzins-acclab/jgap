#ifndef JGAP_SQUAREDEXPKERNEL_HPP
#define JGAP_SQUAREDEXPKERNEL_HPP

#include "Kernel.hpp"
#include <cmath>

namespace jgap {

    template<size_t ExpDimensions, size_t CutoffDimensions>
    class SquaredExpKernel : public Kernel<ExpDimensions + CutoffDimensions> {
    public:
        static constexpr size_t ExpDim = ExpDimensions;
        static constexpr size_t CutoffDim = CutoffDimensions;
        static constexpr size_t TotalDimensions = ExpDimensions + CutoffDimensions;

        using KernelValueAndGradient = Kernel<TotalDimensions>::KernelValueAndGradient;

        SquaredExpKernel() = default;

        SquaredExpKernel(const Real energy_scale, const std::array<Real, ExpDimensions>& length_scales) {
            prefactor = energy_scale * energy_scale;
            for (size_t dim = 0; dim < ExpDimensions; dim++) {
                inverse_length_scales_squared[dim] = 1.0 / (length_scales[dim] * length_scales[dim]);
            }
        }

        // energy_scale and length_scales as passed to the constructor (prefactor = energy_scale^2,
        // inverse_length_scales_squared[d] = 1 / length_scales[d]^2).
        Real getEnergyScale() const {
            return std::sqrt(prefactor);
        }

        std::array<Real, ExpDimensions> getLengthScales() const {
            std::array<Real, ExpDimensions> length_scales{};
            for (size_t dim = 0; dim < ExpDimensions; dim++) {
                length_scales[dim] = 1.0 / std::sqrt(inverse_length_scales_squared[dim]);
            }
            return length_scales;
        }

        Real value(const Descriptor<TotalDimensions> &q1,
                   const Descriptor<TotalDimensions> &q2) const override {

            Real exp_argument = 0;
            for (size_t dim = 0; dim < ExpDimensions; dim++) {
                Real diff = q1.value[dim] - q2.value[dim];
                exp_argument += diff * diff * inverse_length_scales_squared[dim];
            }
            Real val = prefactor * std::exp(-0.5 * exp_argument);

            for (size_t dim = ExpDimensions; dim < TotalDimensions; dim++) {
                val *= q1.value[dim] * q2.value[dim];
            }

            return val;
        }

        KernelValueAndGradient valueAndGradient(const Descriptor<TotalDimensions> &sparse_point,
                                                const Descriptor<TotalDimensions> &q) const override {

            Real exp_argument = 0;
            for (size_t dim = 0; dim < ExpDimensions; dim++) {
                Real diff = q.value[dim] - sparse_point.value[dim];
                exp_argument += diff * diff * inverse_length_scales_squared[dim];
            }
            Real val = prefactor * std::exp(-0.5 * exp_argument);

            for (size_t dim = ExpDimensions; dim < TotalDimensions; dim++) {
                val *= sparse_point.value[dim];
            }

            std::array<Real, TotalDimensions> gradient{};
            for (size_t dim = ExpDimensions; dim < TotalDimensions; dim++) {
                gradient[dim] = val;
            }

            for (size_t dim = ExpDimensions; dim < TotalDimensions; dim++) {
                val *= q.value[dim];

                for (size_t j = ExpDimensions; j < dim; j++) {
                    gradient[j] *= q.value[dim];
                }

                for (size_t j = dim + 1; j < TotalDimensions; j++) {
                    gradient[j] *= q.value[dim];
                }
            }

            for (size_t dim = 0; dim < ExpDimensions; dim++) {
                gradient[dim] = val * (sparse_point.value[dim] - q.value[dim]) * inverse_length_scales_squared[dim];
            }

            return {
                .value = val,
                .gradient = gradient
            };
        }

    private:
        Real prefactor{};
        std::array<Real, ExpDimensions> inverse_length_scales_squared{};
    };
}

#endif