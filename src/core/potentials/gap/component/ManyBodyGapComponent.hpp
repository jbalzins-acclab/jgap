#ifndef JGAP_MANYBODYGAPCOMPONENT_HPP
#define JGAP_MANYBODYGAPCOMPONENT_HPP

#include "core/atomic/energy/AtomicQuantities.hpp"
#include "core/kernels/Kernel.hpp"

#include <memory>
#include "GapComponent.hpp"
#include "core/Matrix.hpp"
#include "core/sparsification/Sparsifier.hpp"
#include "core/transform/manybody/NBodyAggregator.hpp"

namespace jgap {

    template<size_t Dim, CKernelOfDim<Dim> TKernel>
    class ManyBodyGapComponent final : public GapComponent {
    public:
        ManyBodyGapComponent(const ValuePtr<NBodyAggregator<Dim>>& aggregator, const TKernel& kernel,
                             std::vector<Descriptor<Dim>> sparse_points,
                             const std::vector<Real>& optional_coeffs = {}) :
            aggregator(aggregator), kernel(kernel), sparse_points(std::move(sparse_points)) {
            if (!optional_coeffs.empty()) {
                this->setCoefficients(optional_coeffs);
            }
        }

        ManyBodyGapComponent(const ValuePtr<NBodyAggregator<Dim>>& aggregator, const TKernel& kernel,
                             const Sparsifier<Dim>& sparsifier, const std::vector<Atoms>& training_data,
                             const std::vector<Real>& optional_coeffs = {}) :
            ManyBodyGapComponent(aggregator, kernel,
                                 sparsifier.selectSparsePoints(getAllDescriptors(training_data, aggregator)),
                                 optional_coeffs) {}

        std::optional<AtomicQuantities> covariate(const NeighbourLists& nl) const override {
            auto aggregated_descriptors = aggregator->aggregate(nl);
            if (aggregated_descriptors.values.empty()) {
                return std::nullopt;
            }

            AtomicQuantities result(nSparsePoints(), nl.nAtoms());

            for (size_t i = 0; i < aggregated_descriptors.values.size(); i++) {
                const auto& descriptor_value = aggregated_descriptors.values[i];

                for (size_t sparse_idx = 0; sparse_idx < nSparsePoints(); sparse_idx++) {
                    const auto& sparse_desc = sparse_points[sparse_idx];
                    const auto [K, gradK_wrt_q] = kernel.valueAndGradient(sparse_desc, descriptor_value);

                    result.energy(sparse_idx) += K;

                    // Final chain rule combining gradK_wrt_q with the accumulated per-dimension forces and virials
                    for (size_t dim = 0; dim < Dim; dim++) {
                        result.virials(sparse_idx) += aggregated_descriptors.virials[i][dim] * gradK_wrt_q[dim];
                    }

                    for (size_t j = 0; j < nl.nAtoms(); j++) {
                        for (size_t dim = 0; dim < Dim; dim++) {
                            result.force(sparse_idx, j) += aggregated_descriptors.force(i, j)[dim] * gradK_wrt_q[dim];
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
                    result(i, j) = kernel.value(sparse_points[i], sparse_points[j]);
                    result(j, i) = result(i, j);
                }
            }
            return result;
        }

        size_t nSparsePoints() const override { return sparse_points.size(); }

        Cutoffs getCutoffs() const override { return aggregator->getCutoffs(); }

        ManyBodyGapComponent* clone() const override { return new ManyBodyGapComponent(*this); }

        const ValuePtr<NBodyAggregator<Dim>>& getAggregator() const { return aggregator; }
        const TKernel& getKernel() const { return kernel; }
        const std::vector<Descriptor<Dim>>& getSparsePoints() const { return sparse_points; }

        void tabulate(TabulationData& tables) const override {
            const auto& coeffs = this->getCoefficients();
            if (coeffs.empty()) {
                JGAP_LOG_AND_THROW("Coefficients must be set before tabulation");
            }

            if constexpr (Dim == 1) {
                aggregator->tabulateNewManyBodyGrid(tables);
                auto& eam_grids = tables.eam_grids_vec.back();

                for (auto cell: eam_grids.value_grid) {
                    for (size_t sparse_idx = 0; sparse_idx < sparse_points.size(); sparse_idx++) {
                        cell.value += coeffs[sparse_idx] * kernel.value(sparse_points[sparse_idx], cell.pos);
                    }
                }

            } else {
                JGAP_LOG_AND_THROW("Tabulation not implemented");
            }
        }

    private:
        static std::vector<Descriptor<Dim>> getAllDescriptors(const std::vector<Atoms>& training_data,
                                                              const ValuePtr<NBodyAggregator<Dim>>& aggregator) {
            std::vector<Descriptor<Dim>> all_descriptors;
            Real cutoff = aggregator->getCutoffs().maxOverall();
            for (const auto& atoms: training_data) {
                NeighbourLists nl(atoms, cutoff);
                auto aggregated = aggregator->aggregate(nl);
                for (const auto& desc: aggregated.values) {
                    all_descriptors.push_back({desc});
                }
            }
            return all_descriptors;
        }

        ValuePtr<NBodyAggregator<Dim>> aggregator;
        TKernel kernel;

        std::vector<Descriptor<Dim>> sparse_points;
    };
}

#endif
