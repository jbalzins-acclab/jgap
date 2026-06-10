#ifndef JGAP_MANYBODYGAPCOMPONENT_HPP
#define JGAP_MANYBODYGAPCOMPONENT_HPP

#include "core/atomic/energy/AtomicQuantities.hpp"
#include "core/kernels/Kernel.hpp"
#include "../../../transform/aggregated/TransformationAggregator.hpp"
#include "GapComponent.hpp"
#include "core/Matrix.hpp"
#include <memory>

namespace jgap {

    template<size_t Dim>
    class ManyBodyGapComponent final : public GapComponent {
    public:
        ManyBodyGapComponent(std::unique_ptr<TransformationAggregator<Dim>> aggregator,
                             std::unique_ptr<Kernel<Dim>> kernel,
                             std::vector<Descriptor<Dim>> sparse_points,
                             const std::vector<Real>& optional_coeffs = {})
            : aggregator(std::move(aggregator)),
              kernel(std::move(kernel)),
              sparse_points(std::move(sparse_points)) {

            if (!optional_coeffs.empty()) {
                setCoefficients(optional_coeffs);
            }
        }

        std::optional<AtomicQuantities> covariate(const NeighbourList& nl) const override {
            auto aggregated_descriptors = aggregator->aggregate(nl);

            if (aggregated_descriptors.empty()) {
                return std::nullopt;
            }

            AtomicQuantities result(nSparsePoints(), nl.nAtoms());

            for (const auto& [atom_idx, descriptor] : aggregated_descriptors) {
                for (size_t sparse_idx = 0; sparse_idx < nSparsePoints(); sparse_idx++) {
                    const auto& sparse_desc = sparse_points[sparse_idx];
                    const auto [K, gradK_wrt_q] = kernel->valueAndGradient(sparse_desc.value, descriptor.value);

                    result.energy(sparse_idx) += K;

                    // Final chain rule combining gradK_wrt_q with the accumulated per-dimension forces and virials
                    for (size_t dim = 0; dim < Dim; dim++) {
                        result.virials(sparse_idx) += descriptor.virials[dim] * gradK_wrt_q[dim];

                        for (size_t j = 0; j < nl.nAtoms(); j++) {
                            result.force(sparse_idx, j) += descriptor.forces[j][dim] * gradK_wrt_q[dim];
                        }
                    }
                }
            }

            return result;
        }

        Matrix sparseToSparseCovariance() const override {
            Matrix result(nSparsePoints(), nSparsePoints());
            for (size_t i = 0; i < nSparsePoints(); i++) {
                for (size_t j = i; j < nSparsePoints(); j++) {
                    result(i, j) = kernel->value(sparse_points[i].value, sparse_points[j].value);
                    result(j, i) = result(i, j);
                }
            }
            return result;
        }

        size_t nSparsePoints() const override {
            return sparse_points.size();
        }

        Cutoffs getCutoffs() const override {
            return aggregator->getCutoffs();
        }

    public:
        std::unique_ptr<TransformationAggregator<Dim>> aggregator;
        std::unique_ptr<Kernel<Dim>> kernel;
        std::vector<Descriptor<Dim>> sparse_points;
    };
}

#endif