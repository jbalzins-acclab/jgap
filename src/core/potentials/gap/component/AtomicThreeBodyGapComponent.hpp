#ifndef JGAP_ATOMICTHREEBODYGAPCOMPONENT_HPP
#define JGAP_ATOMICTHREEBODYGAPCOMPONENT_HPP

#include <cassert>
#include <ranges>
#include <set>
#include "../../../transform/nbody/3b/ThreeBodyTransformation.hpp"
#include "GapComponent.hpp"
#include "core/Matrix.hpp"
#include "core/atomic/energy/AtomicQuantities.hpp"
#include "core/atomic/iteration/AtomicCluster3Expansion.hpp"
#include "core/kernels/Kernel.hpp"
#include "core/sparsification/Sparsifier.hpp"

namespace jgap {

    template<size_t Dim, CKernelOfDim<Dim> TKernel>
    class AtomicThreeBodyGapComponent : public GapComponent {
    public:
        static constexpr size_t Dependencies = 3;

        static std::vector<AtomicThreeBodyGapComponent> createComponents(
            const std::vector<Atoms>& training_data, const ValuePtr<ThreeBodyTransformation<Dim>>& transformation,
            const TKernel& kernel, const Sparsifier<Dim>& sparsifier) {
            std::set<Species3AtomicSorted> all_species_sets;
            Real cutoff = transformation->getCutoffs().maxOverall();

            for (const auto& atoms: training_data) {
                NeighbourLists nl(atoms, cutoff);
                auto sets = Species3AtomicSorted::getAll(nl);
                all_species_sets.insert(sets.begin(), sets.end());
            }

            std::vector<AtomicThreeBodyGapComponent> components;
            for (const auto& species_set: all_species_sets) {
                components.emplace_back(species_set, transformation, kernel, sparsifier, training_data);
            }

            return components;
        }

        AtomicThreeBodyGapComponent(const Species3AtomicSorted species,
                                    const ValuePtr<ThreeBodyTransformation<Dim>>& transformation, const TKernel& kernel,
                                    const std::vector<Descriptor<Dim>>& sparse_points,
                                    const std::vector<Real>& optional_coeffs = {}) :
            species(species), transformation(transformation), kernel(kernel), sparse_points(sparse_points) {
            if (!optional_coeffs.empty()) {
                this->setCoefficients(optional_coeffs);
            }
        }

        AtomicThreeBodyGapComponent(const Species3AtomicSorted species,
                                    const ValuePtr<ThreeBodyTransformation<Dim>>& transformation, const TKernel& kernel,
                                    const Sparsifier<Dim>& sparsifier, const std::vector<Atoms>& training_data) :
            AtomicThreeBodyGapComponent(
                species, transformation, kernel,
                sparsifier.selectSparsePoints(getAllDescriptors(training_data, species, transformation))) {}

        std::optional<AtomicQuantities> covariate(const NeighbourLists& neighbour_list) const override {
            AtomicQuantities result(nSparsePoints(), neighbour_list.nAtoms());
            bool found_any = false;

            auto it = neighbour_list.atoms_by_species.find(species.root);
            if (it != neighbour_list.atoms_by_species.end()) {
                for (size_t atom_index: it->second) {
                    if (covariate_atom(atom_index, neighbour_list, result)) {
                        found_any = true;
                    }
                }
            }

            if (!found_any) {
                return std::nullopt;
            }

            if (transformation->isSwapInvariant(1, 2)) {
                constexpr Real iteration_reduction_factor = 2.0;
                result *= iteration_reduction_factor;
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

        Cutoffs getCutoffs() const override { return transformation->getCutoffs(); }

        AtomicThreeBodyGapComponent* clone() const override { return new AtomicThreeBodyGapComponent(*this); }

        const Species3AtomicSorted& getSpecies() const { return species; }
        const ValuePtr<ThreeBodyTransformation<Dim>>& getTransformation() const { return transformation; }
        const TKernel& getKernel() const { return kernel; }
        const std::vector<Descriptor<Dim>>& getSparsePoints() const { return sparse_points; }

        void tabulate(TabulationData& tables) const override {
            const auto& coeffs = this->getCoefficients();
            if (coeffs.empty()) {
                JGAP_LOG_AND_THROW("Coefficients must be set before tabulation");
            }

            auto& table = tables.three_body_grids.getValueGrid(species);

            const Real iteration_reduction_factor = transformation->isSwapInvariant(1, 2) ? 2.0 : 1.0;

            for (size_t i = 0; i < table.data_flat.size(); ++i) {
                auto indices = table.getIndices(i);
                auto pos = table.getCoord(indices);

                auto cluster = TabulationData::gridPosAsCluster3(pos);
                auto transformed = transformation->evaluate(cluster);

                Real& value = table.data_flat[i];
                for (size_t sparse_idx = 0; sparse_idx < sparse_points.size(); sparse_idx++) {
                    value += coeffs[sparse_idx] * kernel.value(sparse_points[sparse_idx], transformed) *
                             iteration_reduction_factor;
                }
            }
        }

    private:
        bool covariate_atom(size_t atom_index, const NeighbourLists& neighbour_list, AtomicQuantities& result) const {
            AtomicCluster3Expansion expansion(species);
            auto expansion_result = expansion.find(atom_index, neighbour_list, CalculationType::WithGradients);
            const auto& clusters = expansion_result.clusters;

            if (clusters.empty()) {
                return false;
            }

            assert(expansion_result.derivatives.has_value());

            const bool needs_permutation = !transformation->isSwapInvariant(1, 2);

            for (const auto& [cluster, cluster_derivs]: std::views::zip(clusters, *expansion_result.derivatives)) {
                auto descriptor = transformation->evaluateAndDifferentiate(cluster);

                for (size_t sparse_idx = 0; sparse_idx < sparse_points.size(); sparse_idx++) {
                    auto [K, gradK_wrt_q] = kernel.valueAndGradient(sparse_points[sparse_idx], descriptor.value);

                    result.energy(sparse_idx) += K;

                    for (size_t i = 0; i < 3; i++) {
                        for (size_t j = i + 1; j < 3; j++) {
                            const auto sep_idx = flattenedIndex(i, j);
                            const auto& separation_derivatives = cluster_derivs.derivativesBetween(i, j);
                            const auto& derivatives_wrt_r_norms = descriptor.derivatives[sep_idx];

                            Real dK_drnorm = 0.0;
                            for (size_t dim = 0; dim < Dim; dim++) {
                                dK_drnorm += derivatives_wrt_r_norms[dim] * gradK_wrt_q[dim];
                            }

                            result.force(sparse_idx, cluster.atom_indexes[i]) +=
                                separation_derivatives.direction * dK_drnorm;
                            result.force(sparse_idx, cluster.atom_indexes[j]) -=
                                separation_derivatives.direction * dK_drnorm;

                            result.virials(sparse_idx) += separation_derivatives.virials * dK_drnorm;
                        }
                    }
                }
                if (needs_permutation) {
                    Cluster3 permuted_cluster = cluster;
                    std::swap(permuted_cluster.atom_indexes[1], permuted_cluster.atom_indexes[2]);
                    std::swap(permuted_cluster.separation_magnitudes[0], permuted_cluster.separation_magnitudes[1]);

                    Cluster3Derivatives permuted_derivs = cluster_derivs;
                    std::swap(permuted_derivs.val[0], permuted_derivs.val[1]);
                    permuted_derivs.val[2].direction *= -1.0;

                    auto descriptor2 = transformation->evaluateAndDifferentiate(permuted_cluster);

                    for (size_t sparse_idx = 0; sparse_idx < sparse_points.size(); sparse_idx++) {
                        auto [K, gradK_wrt_q] = kernel.valueAndGradient(sparse_points[sparse_idx], descriptor2.value);

                        result.energy(sparse_idx) += K;

                        for (size_t i = 0; i < 3; i++) {
                            for (size_t j = i + 1; j < 3; j++) {
                                const auto sep_idx = flattenedIndex(i, j);
                                const auto& separation_derivatives = permuted_derivs.derivativesBetween(i, j);
                                const auto& derivatives_wrt_r_norms = descriptor2.derivatives[sep_idx];

                                Real dK_drnorm = 0.0;
                                for (size_t dim = 0; dim < Dim; dim++) {
                                    dK_drnorm += derivatives_wrt_r_norms[dim] * gradK_wrt_q[dim];
                                }

                                result.force(sparse_idx, permuted_cluster.atom_indexes[i]) +=
                                    separation_derivatives.direction * dK_drnorm;
                                result.force(sparse_idx, permuted_cluster.atom_indexes[j]) -=
                                    separation_derivatives.direction * dK_drnorm;

                                result.virials(sparse_idx) += separation_derivatives.virials * dK_drnorm;
                            }
                        }
                    }
                }
            }
            return true;
        }

        static std::vector<Descriptor<Dim>> getAllDescriptors(
            const std::vector<Atoms>& training_data, const Species3AtomicSorted& species,
            const ValuePtr<ThreeBodyTransformation<Dim>>& transformation) {
            std::vector<Descriptor<Dim>> all_descriptors;
            AtomicCluster3Expansion expansion(species);
            for (const auto& atoms: training_data) {
                NeighbourLists nl(atoms, transformation->getCutoffs().maxOverall());

                auto it = nl.atoms_by_species.find(species.root);
                if (it != nl.atoms_by_species.end()) {
                    for (size_t atom_index: it->second) {
                        auto clusters = expansion.find(atom_index, nl, CalculationType::ValueOnly).clusters;
                        for (const auto& cluster: clusters) {
                            all_descriptors.push_back(transformation->evaluate(cluster));
                        }
                    }
                }
            }
            return all_descriptors;
        }

        Species3AtomicSorted species;
        ValuePtr<ThreeBodyTransformation<Dim>> transformation;
        TKernel kernel;
        std::vector<Descriptor<Dim>> sparse_points;
    };
}

#endif
