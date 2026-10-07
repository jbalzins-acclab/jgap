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

        CauchyKernel(const double energy_scale, const std::array<double, ExpDimensions>& length_scales) {
            prefactor = energy_scale * energy_scale;
            for (size_t dim = 0; dim < ExpDimensions; dim++) {
                inverse_length_scales_squared[dim] = 1.0 / (length_scales[dim] * length_scales[dim]);
            }
        }

        double getEnergyScale() const { return std::sqrt(prefactor); }

        std::array<double, ExpDimensions> getLengthScales() const {
            std::array<double, ExpDimensions> length_scales{};
            for (size_t dim = 0; dim < ExpDimensions; dim++) {
                length_scales[dim] = 1.0 / std::sqrt(inverse_length_scales_squared[dim]);
            }
            return length_scales;
        }

        KernelValueAndGradient valueAndGradient(const Descriptor<TotalDimensions>& sparse_point,
                                                const Descriptor<TotalDimensions>& q) const override {
            double dist_sq = 0.0;
            for (size_t dim = 0; dim < ExpDimensions; dim++) {
                double diff = q[dim] - sparse_point[dim];
                dist_sq += diff * diff * inverse_length_scales_squared[dim];
            }
            double denom = 1.0 + dist_sq;
            double inv_denom = 1.0 / denom;
            double base_val = prefactor * inv_denom;

            double val = base_val;
            std::array<double, TotalDimensions> gradient{};

            if constexpr (CutoffDimensions == 1) {
                gradient[ExpDimensions] = base_val * sparse_point[ExpDimensions];
                val = base_val * sparse_point[ExpDimensions] * q[ExpDimensions];
            }

            double factor = val * inv_denom * 2.0;
            for (size_t dim = 0; dim < ExpDimensions; dim++) {
                gradient[dim] = factor * (sparse_point[dim] - q[dim]) * inverse_length_scales_squared[dim];
            }

            return {.value = val, .gradient = gradient};
        }

        CauchyKernel* clone() const override { return new CauchyKernel(*this); }

    private:
        double prefactor{};
        std::array<double, ExpDimensions> inverse_length_scales_squared{};
    };
}

#endif
