#ifndef JGAP_MANYBODYGAPCOMPONENT_HPP
#define JGAP_MANYBODYGAPCOMPONENT_HPP

#include "core/atomic/energy/AtomicQuantities.hpp"
#include "core/kernels/Kernel.hpp"
#include "../../../transform/aggregated/TransformationAggregator.hpp"
#include "core/sparsification/Sparsifier.hpp"
#include "GapComponent.hpp"
#include "core/RowMajorMatrix.hpp"
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
                this->setCoefficients(optional_coeffs);
            }
        }

        ManyBodyGapComponent(const ValuePtr<TransformationAggregator<Dim>>& aggregator,
                             const TKernel& kernel,
                             const Sparsifier<Dim>& sparsifier,
                             const std::vector<Atoms>& training_data,
                             const std::vector<Real>& optional_coeffs = {})
            : ManyBodyGapComponent(aggregator, kernel,
                                   sparsifier.selectSparsePoints(getAllDescriptors(training_data, aggregator)),
                                   optional_coeffs) {
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

        RowMajorMatrix sparseToSparseCovariance() const override {
            RowMajorMatrix result(nSparsePoints(), nSparsePoints());
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

        ManyBodyGapComponent* clone() const override {
            return new ManyBodyGapComponent(*this);
        }

        const ValuePtr<TransformationAggregator<Dim>>& getAggregator() const { return aggregator; }
        const TKernel& getKernel() const { return kernel; }
        const std::vector<Descriptor<Dim>>& getSparsePoints() const { return sparse_points; }

        void tabulate(TabulationData &tables) const override {
            const auto& coeffs = this->getCoefficients();
            if (coeffs.empty()) {
                JGAP_LOG_AND_THROW("Coefficients must be set before tabulation");
            }

            if constexpr (Dim == 1) {

                aggregator->tabulateNewManyBodyGrid(tables);
                auto& eam_grids = tables.eam_grids_vec.back();

                for (auto cell: eam_grids.value_grid) {
                    for (size_t sparse_idx = 0; sparse_idx < sparse_points.size(); ++sparse_idx) {
                        cell.value += coeffs[sparse_idx] * kernel.value(sparse_points[sparse_idx].value, cell.pos);
                    }
                }

            } else {
                JGAP_LOG_AND_THROW("Tabulation not implemented");
            }
        }

    private:
        static std::vector<Descriptor<Dim>> getAllDescriptors(const std::vector<Atoms>& training_data,
                                                              const ValuePtr<TransformationAggregator<Dim>>& aggregator) {
            std::vector<Descriptor<Dim>> all_descriptors;
            Real cutoff = aggregator->getCutoffs().maxOverall();
            for (const auto& atoms : training_data) {
                NeighbourList nl(atoms, cutoff);
                auto aggregated = aggregator->aggregate(nl);
                for (const auto& [idx, desc] : aggregated) {
                    all_descriptors.push_back({desc.value});
                }
            }
            return all_descriptors;
        }

        ValuePtr<TransformationAggregator<Dim>> aggregator;
        TKernel kernel;

        std::vector<Descriptor<Dim>> sparse_points;
    };
}

#endif