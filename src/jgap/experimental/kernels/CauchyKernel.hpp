#ifndef JGAP_CAUCHYKERNEL_HPP
#define JGAP_CAUCHYKERNEL_HPP

#include <cmath>
#include "jgap/core/kernels/Kernel.hpp"

namespace jgap {

    template<size_t ExpDimensions, size_t CutoffDimensions>
        requires(CutoffDimensions <= 1)
    class CauchyKernel : public Kernel<ExpDimensions + CutoffDimensions> {
    public:
        static constexpr size_t ExpDim = ExpDimensions;
        static constexpr size_t CutoffDim = CutoffDimensions;
        static constexpr size_t TotalDimensions = ExpDimensions + CutoffDimensions;

        using KernelValueAndGradient = Kernel<TotalDimensions>::KernelValueAndGradient;

        CauchyKernel() = default;

        CauchyKernel(const Real energy_scale, const std::array<Real, ExpDimensions>& length_scales) {
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

        KernelValueAndGradient valueAndGradient(const Descriptor<TotalDimensions>& sparse_point,
                                                const Descriptor<TotalDimensions>& q) const override {
            Real dist_sq = 0.0_r;
            for (size_t dim = 0; dim < ExpDimensions; dim++) {
                Real diff = q[dim] - sparse_point[dim];
                dist_sq += diff * diff * inverse_length_scales_squared[dim];
            }
            Real denom = 1.0_r + dist_sq;
            Real inv_denom = 1.0_r / denom;
            Real base_val = prefactor * inv_denom;

            Real val = base_val;
            std::array<Real, TotalDimensions> gradient{};

            if constexpr (CutoffDimensions == 1) {
                gradient[ExpDimensions] = base_val * sparse_point[ExpDimensions];
                val = base_val * sparse_point[ExpDimensions] * q[ExpDimensions];
            }

            Real factor = val * inv_denom * 2.0_r;
            for (size_t dim = 0; dim < ExpDimensions; dim++) {
                gradient[dim] = factor * (sparse_point[dim] - q[dim]) * inverse_length_scales_squared[dim];
            }

            return {.value = val, .gradient = gradient};
        }

        CauchyKernel* clone() const override { return new CauchyKernel(*this); }

    private:
        Real prefactor{};
        std::array<Real, ExpDimensions> inverse_length_scales_squared{};
    };
}

#endif
