#ifndef JGAP_MANYBODYGAPCOMPONENT_HPP
#define JGAP_MANYBODYGAPCOMPONENT_HPP

#include "core/atomic/energy/AtomicQuantities.hpp"
#include "core/kernels/Kernel.hpp"
#include "../../../transform/aggregated/TransformationAggregator.hpp"
#include "GapComponent.hpp"
#include "core/Matrix.hpp"
#include <memory>

namespace jgap {

    template<size_t Dim, CKernelOfDim<Dim> TKernel>
    class ManyBodyGapComponent final : public GapComponent {
    public:
        ManyBodyGapComponent(const ValuePtr<TransformationAggregator<Dim>>& aggregator,
                             const TKernel& kernel,
                             std::vector<Descriptor<Dim>> sparse_points,
                             const std::vector<Real>& optional_coeffs = {})
            : aggregator(aggregator),
              kernel(kernel),
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
                    const auto [K, gradK_wrt_q] = kernel.valueAndGradient(sparse_desc.value, descriptor.value);

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
                    result(i, j) = kernel.value(sparse_points[i].value, sparse_points[j].value);
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

        std::unique_ptr<GapComponent> clone() const override {
            return std::make_unique<ManyBodyGapComponent>(*this);
        }

        void tabulate(TabulationData &tables) const override {
            if constexpr (Dim == 1) {

                aggregator->tabulateNewManyBodyGrid(tables);
                auto& eam_grids = tables.eam_grids_vec.back();

                for (auto cell: eam_grids.value_grid) {
                    for (const auto& sparse_point: sparse_points) {
                        cell.value += kernel.value(sparse_point.value, cell.pos);
                    }
                }

            } else {
                JGAP_LOG_AND_THROW("Tabulation not implemented");
            }
        }

    private:
        ValuePtr<TransformationAggregator<Dim>> aggregator;
        TKernel kernel;

        std::vector<Descriptor<Dim>> sparse_points;
    };
}

#endif