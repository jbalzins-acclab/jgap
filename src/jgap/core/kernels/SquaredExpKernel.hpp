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

        SquaredExpKernel(const double energy_scale, const std::array<double, ExpDimensions>& length_scales) {
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

        double value(const Descriptor<TotalDimensions>& q1, const Descriptor<TotalDimensions>& q2) const override {
            return Kernel<TotalDimensions>::value(q1, q2);
        }

        KernelValueAndGradient valueAndGradient(
            const Descriptor<TotalDimensions>& sparse_point, const Descriptor<TotalDimensions>& q
        ) const override {
            double exp_argument = 0.0;
            for (size_t dim = 0; dim < ExpDimensions; dim++) {
                double diff = q[dim] - sparse_point[dim];
                exp_argument += diff * diff * inverse_length_scales_squared[dim];
            }
            double val = prefactor * std::exp(-0.5 * exp_argument);

            std::array<double, TotalDimensions> gradient{};

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
        double prefactor{};
        std::array<double, ExpDimensions> inverse_length_scales_squared{};
    };
}

#endif
