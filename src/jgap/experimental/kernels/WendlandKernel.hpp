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

        WendlandKernel(const Real energy_scale, const std::array<Real, ExpDimensions>& length_scales) {
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
            Real dist_sq = 0.0_r;
            for (size_t dim = 0; dim < ExpDimensions; dim++) {
                Real diff = q1[dim] - q2[dim];
                dist_sq += diff * diff * inverse_length_scales_squared[dim];
            }

            if (dist_sq >= 1.0_r) {
                return 0.0_r;
            }

            Real r = std::sqrt(dist_sq);
            Real omr = 1.0_r - r;
            Real omr2 = omr * omr;
            Real omr4 = omr2 * omr2;
            Real val = prefactor * omr4 * (1.0_r + 4.0_r * r);

            if constexpr (CutoffDimensions == 1) {
                val *= q1[ExpDimensions] * q2[ExpDimensions];
            }

            return val;
        }

        KernelValueAndGradient valueAndGradient(
            const Descriptor<TotalDimensions>& sparse_point, const Descriptor<TotalDimensions>& q
        ) const override {
            Real dist_sq = 0.0_r;
            for (size_t dim = 0; dim < ExpDimensions; dim++) {
                Real diff = q[dim] - sparse_point[dim];
                dist_sq += diff * diff * inverse_length_scales_squared[dim];
            }

            Real val = 0.0_r;
            std::array<Real, TotalDimensions> gradient{};

            if (dist_sq < 1.0_r) {
                Real r = std::sqrt(dist_sq);
                Real omr = 1.0_r - r;
                Real omr2 = omr * omr;
                Real omr3 = omr2 * omr;
                Real omr4 = omr2 * omr2;

                Real base_val = prefactor * omr4 * (1.0_r + 4.0_r * r);
                val = base_val;

                if constexpr (CutoffDimensions == 1) {
                    gradient[ExpDimensions] = base_val * sparse_point[ExpDimensions];
                    val = base_val * sparse_point[ExpDimensions] * q[ExpDimensions];
                }

                Real factor = 20.0_r * prefactor * omr3;
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
        Real prefactor{};
        std::array<Real, ExpDimensions> inverse_length_scales_squared{};
    };
}

#endif
