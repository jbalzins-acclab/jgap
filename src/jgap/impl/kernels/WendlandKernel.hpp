#ifndef JGAP_WENDLANDKERNEL_HPP
#define JGAP_WENDLANDKERNEL_HPP

#include <algorithm>
#include <cmath>
#include "jgap/core/kernels/Kernel.hpp"

namespace jgap {

    template<size_t ExpDimensions, size_t CutoffDimensions>
        requires(CutoffDimensions <= 1)
    class WendlandKernel : public Kernel<ExpDimensions + CutoffDimensions> {
    public:
        static constexpr size_t ExpDim = ExpDimensions;
        static constexpr size_t CutoffDim = CutoffDimensions;
        static constexpr size_t TotalDimensions = ExpDimensions + CutoffDimensions;

        using KernelValueAndGradient = Kernel<TotalDimensions>::KernelValueAndGradient;

        WendlandKernel() = default;

        WendlandKernel(const double energy_scale, const std::array<double, ExpDimensions>& length_scales) {
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

        KernelValueAndGradient valueAndGradient(
            const Descriptor<TotalDimensions>& sparse_point, const Descriptor<TotalDimensions>& q
        ) const override {
            double dist_sq = 0.0;
            for (size_t dim = 0; dim < ExpDimensions; dim++) {
                double diff = q[dim] - sparse_point[dim];
                dist_sq += diff * diff * inverse_length_scales_squared[dim];
            }

            double val = 0.0;
            std::array<double, TotalDimensions> gradient{};

            if (dist_sq < 1.0) {
                double r = std::sqrt(dist_sq);
                double omr = 1.0 - r;
                double omr2 = omr * omr;
                double omr3 = omr2 * omr;
                double omr4 = omr2 * omr2;

                double base_val = prefactor * omr4 * (1.0 + 4.0 * r);
                val = base_val;

                if constexpr (CutoffDimensions == 1) {
                    gradient[ExpDimensions] = base_val * sparse_point[ExpDimensions];
                    val = base_val * sparse_point[ExpDimensions] * q[ExpDimensions];
                }

                double factor = 20.0 * prefactor * omr3;
                if constexpr (CutoffDimensions == 1) {
                    factor *= sparse_point[ExpDimensions] * q[ExpDimensions];
                }

                for (size_t dim = 0; dim < ExpDimensions; dim++) {
                    gradient[dim] = factor * (sparse_point[dim] - q[dim]) * inverse_length_scales_squared[dim];
                }
            }

            return {.value = val, .gradient = gradient};
        }

        WendlandKernel* clone() const override { return new WendlandKernel(*this); }

    private:
        double prefactor{};
        std::array<double, ExpDimensions> inverse_length_scales_squared{};
    };
}

#endif
