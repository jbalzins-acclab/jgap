#ifndef JGAP_TWOBODYGAPCOMPONENT_HPP
#define JGAP_TWOBODYGAPCOMPONENT_HPP

#include <cassert>
#include <ranges>
#include <set>
#include "GapComponent.hpp"
#include "core/Matrix.hpp"
#include "core/ValuePtr.hpp"
#include "core/atomic/energy/AtomicQuantities.hpp"
#include "core/atomic/iteration/Cluster2Expansion.hpp"
#include "core/atomic/neighbours/NeighbourLists.hpp"
#include "core/kernels/Kernel.hpp"
#include "core/sparsification/Sparsifier.hpp"
#include "core/transform/nbody/2b/TwoBodyTransformation.hpp"

namespace jgap {

    template<size_t Dim, CKernelOfDim<Dim> TKernel>
    class TwoBodyGapComponent : public GapComponent {
    public:
        static constexpr size_t Dependencies = 1;

        static std::vector<TwoBodyGapComponent> createComponents(
            const std::vector<Atoms>& training_data, const ValuePtr<TwoBodyTransformation<Dim>>& transformation,
            const TKernel& kernel, const Sparsifier<Dim>& sparsifier) {
            std::set<Species2Sorted> all_species_sets;
            Real cutoff = transformation->getCutoffs().maxOverall();

            for (const auto& atoms: training_data) {
                NeighbourLists nl(atoms, cutoff);
                auto sets = Species2Sorted::getAll(nl);
                all_species_sets.insert(sets.begin(), sets.end());
            }

            std::vector<TwoBodyGapComponent> components;
            for (const auto& species_set: all_species_sets) {
                components.emplace_back(species_set, transformation, kernel, sparsifier, training_data);
            }

            return components;
        }

        TwoBodyGapComponent(const Species2Sorted& species, ValuePtr<TwoBodyTransformation<Dim>> transformation,
                            TKernel kernel, std::vector<Descriptor<Dim>> sparse_points,
                            const std::vector<Real>& optional_coeffs = {}) :
            expansion(species),
            transformation(std::move(transformation)),
            kernel(std::move(kernel)),
            sparse_points(std::move(sparse_points)) {
            if (!optional_coeffs.empty()) {
                this->setCoefficients(optional_coeffs);
            }
        }

        TwoBodyGapComponent(const Species2Sorted species, const ValuePtr<TwoBodyTransformation<Dim>>& transformation,
                            const TKernel& kernel, const Sparsifier<Dim>& sparsifier,
                            const std::vector<Atoms>& training_data) :
            expansion(species), transformation(transformation), kernel(kernel) {
            auto all_descriptors = getAllDescriptors(training_data, species, transformation);
            sparse_points = sparsifier.selectSparsePoints(all_descriptors);
        }

        std::optional<AtomicQuantities> covariate(const NeighbourLists& neighbour_list) const override {
            auto expansion_result = expansion.find(neighbour_list, CalculationType::WithGradients);
            assert(expansion_result.derivatives.has_value());
            const auto& clusters = expansion_result.clusters;

            if (clusters.empty()) {
                return std::nullopt;
            }

            AtomicQuantities result(nSparsePoints(), neighbour_list.nAtoms());

            for (const auto& [cluster, cluster_derivs]: std::views::zip(clusters, *expansion_result.derivatives)) {
                auto descriptor = transformation->evaluateAndDifferentiate(cluster);

                for (size_t sparse_idx = 0; sparse_idx < sparse_points.size(); sparse_idx++) {
                    auto [K, gradK_wrt_q] = kernel.valueAndGradient(sparse_points[sparse_idx], descriptor.value);

                    result.energy(sparse_idx) += K;

                    const auto& separation_derivatives = cluster_derivs.dr01;
                    const auto& derivatives_wrt_r_norms = descriptor.derivatives;

                    Real dK_drnorm = 0.0;
                    for (size_t dim = 0; dim < Dim; dim++) {
                        dK_drnorm += derivatives_wrt_r_norms[dim] * gradK_wrt_q[dim];
                    }

                    result.force(sparse_idx, cluster.atom_indexes[0]) += separation_derivatives.direction * dK_drnorm;
                    result.force(sparse_idx, cluster.atom_indexes[1]) -= separation_derivatives.direction * dK_drnorm;

                    result.virials(sparse_idx) += separation_derivatives.virials * dK_drnorm;
                }
            }

            constexpr Real iteration_reduction_factor = 2.0;
            result *= iteration_reduction_factor;

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

        Cutoffs getCutoffs() const override { return transformation->getCutoffs(); }

        TwoBodyGapComponent* clone() const override { return new TwoBodyGapComponent(*this); }

        const Species2Sorted& getSpecies() const { return expansion.getSpecies(); }
        const ValuePtr<TwoBodyTransformation<Dim>>& getTransformation() const { return transformation; }
        const TKernel& getKernel() const { return kernel; }
        const std::vector<Descriptor<Dim>>& getSparsePoints() const { return sparse_points; }

        void tabulate(TabulationData& tables) const override {
            const auto& coeffs = this->getCoefficients();
            if (coeffs.empty()) {
                JGAP_LOG_AND_THROW("Coefficients must be set before tabulation");
            }

            Grid<Dependencies>* table_ref = &tables.two_body_grids.getValueGrid(expansion.getSpecies());
            auto& table = *table_ref;

            for (size_t i = 0; i < table.data_flat.size(); ++i) {
                auto indices = table.getIndices(i);
                auto pos = table.getCoord(indices);

                auto cluster = TabulationData::gridPosAsCluster2(pos);
                auto transformed = transformation->evaluate(cluster);

                Real& value = table.data_flat[i];
                for (size_t sparse_idx = 0; sparse_idx < sparse_points.size(); sparse_idx++) {
                    constexpr Real iteration_reduction_factor = 2.0;
                    value += coeffs[sparse_idx] * kernel.value(sparse_points[sparse_idx], transformed) *
                             iteration_reduction_factor;
                }
            }
        }

    private:
        static std::vector<Descriptor<Dim>> getAllDescriptors(
            const std::vector<Atoms>& training_data, const Species2Sorted& species,
            const ValuePtr<TwoBodyTransformation<Dim>>& transformation) {
            std::vector<Descriptor<Dim>> all_descriptors;
            Cluster2Expansion expansion(species);

            for (const auto& atoms: training_data) {
                NeighbourLists nl(atoms, transformation->getCutoffs().maxOverall());
                auto clusters = expansion.find(nl, CalculationType::ValueOnly).clusters;
                for (const auto& cluster: clusters) {
                    all_descriptors.push_back(transformation->evaluate(cluster));
                }
            }
            return all_descriptors;
        }

        Cluster2Expansion expansion;
        ValuePtr<TwoBodyTransformation<Dim>> transformation;
        TKernel kernel;
        std::vector<Descriptor<Dim>> sparse_points;
    };
}

#endif
