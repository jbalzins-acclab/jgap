#ifndef JGAP_SQUAREDEXPKERNEL_HPP
#define JGAP_SQUAREDEXPKERNEL_HPP

#include <cmath>
#include "Kernel.hpp"

namespace jgap {

    template<size_t ExpDimensions, size_t CutoffDimensions>
        requires(CutoffDimensions <= 1)
    class SquaredExpKernel final : public Kernel<ExpDimensions + CutoffDimensions> {
    public:
        static constexpr size_t ExpDim = ExpDimensions;
        static constexpr size_t CutoffDim = CutoffDimensions;
        static constexpr size_t TotalDimensions = ExpDimensions + CutoffDimensions;

        using KernelValueAndGradient = Kernel<TotalDimensions>::KernelValueAndGradient;

        SquaredExpKernel() = default;

        SquaredExpKernel(const Real energy_scale, const std::array<Real, ExpDimensions>& length_scales) {
            prefactor = energy_scale * energy_scale;
            for (size_t dim = 0; dim < ExpDimensions; dim++) {
                inverse_length_scales_squared[dim] = 1.0_r / (length_scales[dim] * length_scales[dim]);
            }
        }

        Real getEnergyScale() const { return std::sqrt(prefactor); }

        std::array<Real, ExpDimensions> getLengthScales() const {
            std::array<Real, ExpDimensions> length_scales{};
            for (size_t dim = 0; dim < ExpDimensions; dim++) {
                length_scales[dim] = 1.0_r / std::sqrt(inverse_length_scales_squared[dim]);
            }
            return length_scales;
        }

        Real value(const Descriptor<TotalDimensions>& q1, const Descriptor<TotalDimensions>& q2) const override {
            return Kernel<TotalDimensions>::value(q1, q2);
        }

        KernelValueAndGradient valueAndGradient(
            const Descriptor<TotalDimensions>& sparse_point, const Descriptor<TotalDimensions>& q
        ) const override {
            Real exp_argument = 0.0_r;
            for (size_t dim = 0; dim < ExpDimensions; dim++) {
                Real diff = q[dim] - sparse_point[dim];
                exp_argument += diff * diff * inverse_length_scales_squared[dim];
            }
            Real val = prefactor * std::exp(-0.5_r * exp_argument);

            std::array<Real, TotalDimensions> gradient{};

            if constexpr (CutoffDimensions == 1) {
                gradient[ExpDimensions] = val * sparse_point[ExpDimensions];
                val *= sparse_point[ExpDimensions] * q[ExpDimensions];
            }

            for (size_t dim = 0; dim < ExpDimensions; dim++) {
                gradient[dim] = val * (sparse_point[dim] - q[dim]) * inverse_length_scales_squared[dim];
            }

            return {.value = val, .gradient = gradient};
        }

        SquaredExpKernel* clone() const override { return new SquaredExpKernel(*this); }

    private:
        Real prefactor{};
        std::array<Real, ExpDimensions> inverse_length_scales_squared{};
    };
}

#endif
