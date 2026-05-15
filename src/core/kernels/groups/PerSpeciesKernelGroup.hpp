#ifndef JGAP_PERSPECIESKERNELGROUP_HPP
#define JGAP_PERSPECIESKERNELGROUP_HPP

#include "core/atomic/energy/AtomicQuantity.hpp"
#include "core/kernels/Kernel.hpp"
#include "core/transform/Transformer.hpp"
#include "core/sparsification/Sparsifier.hpp"
#include "io/log/CurrentLogger.hpp"
#include <concepts>
#include <optional>

#include "KernelGroup.hpp"
#include "core/Matrix.hpp"

namespace jgap {

    template<typename TKernel>
    requires std::derived_from<TKernel, Kernel<TKernel::Dim, TKernel::Dependencies>>
    class PerSpeciesKernelGroup : public KernelGroup<TKernel::Dim, TKernel::Dependencies> {
    public:
        constexpr static size_t DescriptorDim = TKernel::Dim;
        constexpr static size_t DescriptorDependencies = TKernel::Dependencies;

        using TDescriptor = Descriptor<DescriptorDim, DescriptorDependencies, CalculationType::WithGradients>;

        using TSparsePoint = std::array<Real, DescriptorDim>;
        using THyperParamGroup = std::pair<TKernel, std::vector<TSparsePoint>>;

        PerSpeciesKernelGroup(std::map<SpeciesSet, THyperParamGroup> species_groups);

        static PerSpeciesKernelGroup constructFromTrainingData(
            const std::shared_ptr<Transformer<DescriptorDim, DescriptorDependencies>>& transformer,
            const std::vector<NeighbourList>& training_data,
            const std::map<SpeciesSet, std::shared_ptr<Sparsifier<DescriptorDim>>>& sparsifiers = {},
            std::shared_ptr<Sparsifier<DescriptorDim>> default_sparsifier = nullptr,
            const std::map<SpeciesSet, TKernel>& kernels = {},
            std::optional<TKernel> default_kernel = std::nullopt,
            bool ignore_missing = false);

        std::vector<AtomicQuantity> covariate(
            size_t n_atoms, const std::map<SpeciesSet, std::vector<TDescriptor>> &descriptors) const override;
        std::vector<Matrix> sparseToSparseCovariance() const override;

        size_t nSparsePoints() const override {
            return n_sparse_points;
        }

    private:
        size_t n_sparse_points;
        std::map<SpeciesSet, THyperParamGroup> species_groups;
        std::map<SpeciesSet, size_t> species_group_sizes;
    };

    template<typename TKernel>
    requires std::derived_from<TKernel, Kernel<TKernel::Dim, TKernel::Dependencies>>
    PerSpeciesKernelGroup<TKernel>::PerSpeciesKernelGroup(
        std::map<SpeciesSet, THyperParamGroup> species_groups): species_groups(std::move(species_groups)) {

        n_sparse_points = 0;

        for (const auto& [species, hp_group]: this->species_groups) {
            species_group_sizes[species] = hp_group.second.size();
            n_sparse_points += species_group_sizes[species];
        }
    }

    template<typename TKernel>
    requires std::derived_from<TKernel, Kernel<TKernel::Dim, TKernel::Dependencies>>
    PerSpeciesKernelGroup<TKernel> PerSpeciesKernelGroup<TKernel>::constructFromTrainingData(
        const std::shared_ptr<Transformer<DescriptorDim, DescriptorDependencies>>& transformer,
        const std::vector<NeighbourList>& training_data,
        const std::map<SpeciesSet, std::shared_ptr<Sparsifier<DescriptorDim>>>& sparsifiers,
        std::shared_ptr<Sparsifier<DescriptorDim>> default_sparsifier,
        const std::map<SpeciesSet, TKernel>& kernels,
        std::optional<TKernel> default_kernel,
        bool ignore_missing) {

        std::map<SpeciesSet, std::vector<Descriptor<DescriptorDim>>> all_descriptors;

        for (const auto& nl : training_data) {
            auto batch_descriptors = transformer->transform(nl);
            for (const auto& [species, descriptors] : batch_descriptors) {
                for (const auto& desc : descriptors) {
                    Descriptor<DescriptorDim> flat_desc;
                    flat_desc.value = desc.value;
                    all_descriptors[species].push_back(flat_desc);
                }
            }
        }

        std::map<SpeciesSet, THyperParamGroup> species_groups;

        for (const auto& [species, descriptors] : all_descriptors) {

            std::shared_ptr<Sparsifier<DescriptorDim>> sparsifier_to_use = nullptr;
            if (sparsifiers.contains(species)) {
                sparsifier_to_use = sparsifiers.at(species);
            } else if (default_sparsifier != nullptr) {
                sparsifier_to_use = default_sparsifier;
            } else if (ignore_missing) {
                continue;
            } else {
                JGAP_LOG_AND_THROW("No sparsifier specified for species set '{}' and no default sparsifier provided.", species.toString());
            }

            auto sparse_descriptors = sparsifier_to_use->selectSparsePoints(descriptors);

            std::vector<TSparsePoint> sparse_points;
            sparse_points.reserve(sparse_descriptors.size());
            for (const auto& desc : sparse_descriptors) {
                sparse_points.push_back(desc.value);
            }

            const TKernel* kernel_to_use = nullptr;
            if (kernels.contains(species)) {
                kernel_to_use = &kernels.at(species);
            } else if (default_kernel.has_value()) {
                kernel_to_use = &default_kernel.value();
            } else if (ignore_missing) {
                continue;
            } else {
                JGAP_LOG_AND_THROW("No kernel specified for species set '{}' and no default kernel provided.", species.toString());
            }

            species_groups[species] = { *kernel_to_use, std::move(sparse_points) };
        }

        return PerSpeciesKernelGroup(std::move(species_groups));
    }

    template<typename TKernel>
    requires std::derived_from<TKernel, Kernel<TKernel::Dim, TKernel::Dependencies>>
    std::vector<AtomicQuantity> PerSpeciesKernelGroup<TKernel>::covariate(
        const size_t n_atoms, const std::map<SpeciesSet, std::vector<TDescriptor>> &descriptors) const {

        std::vector results(n_sparse_points, AtomicQuantity());

        size_t results_idx = 0;
        for (const auto& [species, hp_group]: species_groups) {
            if (!descriptors.contains(species)) {
                results_idx += species_group_sizes.at(species);
                continue;
            }

            for (size_t tmp_idx = results_idx; tmp_idx < results_idx + species_group_sizes.at(species); tmp_idx++) {
                results[tmp_idx].forces.resize(n_atoms);
            }

            const auto& species_descriptors = descriptors.at(species);
            const auto& [parametrized_kernel, sparse_points] = hp_group;

            for (const auto& sparse_point: sparse_points) {
                parametrized_kernel.covariance(results[results_idx++], sparse_point, species_descriptors);
            }
        }

        return results;
    }

    template<typename TKernel>
    requires std::derived_from<TKernel, Kernel<TKernel::Dim, TKernel::Dependencies>>
    std::vector<Matrix> PerSpeciesKernelGroup<TKernel>::sparseToSparseCovariance() const {

        std::vector<Matrix> result;

        for (const auto& [species, hp_group]: species_groups) {
            const auto& [parametrized_kernel, sparse_points] = hp_group;

            result.emplace_back(sparse_points.size(), sparse_points.size());
            auto& current_matrix = result.back();
            size_t m = sparse_points.size();

            for (size_t i = 0; i < m; i++) {
                for (size_t j = i; j < m; j++) {
                    current_matrix(i, j) = parametrized_kernel.value(sparse_points[i], sparse_points[j]);
                    current_matrix(j, i) = current_matrix(i, j);
                }
            }
        }

        return result;
    }
}

#endif