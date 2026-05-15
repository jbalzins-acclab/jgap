#ifndef JGAP_GENERICKERNELGROUP_HPP
#define JGAP_GENERICKERNELGROUP_HPP

#include "core/atomic/energy/AtomicQuantity.hpp"
#include "core/kernels/Kernel.hpp"
#include "core/transform/Transformer.hpp"
#include <concepts>
#include <type_traits>

#include "KernelGroup.hpp"
#include "core/Matrix.hpp"

namespace jgap {

    // would be nice to just pass a single pointer to a kernel + hyperparam setter, but then concurrency breaks:(
    template<typename TKernel>
    requires std::derived_from<TKernel, Kernel<TKernel::Dim, TKernel::Dependencies>>
    class GenericKernelGroup : public KernelGroup<TKernel::Dim, TKernel::Dependencies> {
    public:
        constexpr static size_t DescriptorDim = TKernel::Dim;
        constexpr static size_t DescriptorDependencies = TKernel::Dependencies;

        using TDescriptor = Descriptor<DescriptorDim, DescriptorDependencies, CalculationType::WithGradients>;

        using TSparsePoint = std::array<Real, DescriptorDim>;
        using THyperParamGroup = std::pair<TKernel, std::vector<TSparsePoint>>;
        using TSpeciesGroup = std::vector<THyperParamGroup>;

        //GenericKernelGroup() : n_sparse_points(0), species_groups({}), species_group_sizes({}) {}
        GenericKernelGroup(std::map<SpeciesSet, TSpeciesGroup> species_groups);

        std::vector<AtomicQuantity> covariate(size_t n_atoms,
            const std::map<SpeciesSet, std::vector<TDescriptor>> &descriptors) const override;
        std::vector<Matrix> sparseToSparseCovariance() const override;

        size_t nSparsePoints() const override {
            return n_sparse_points;
        }

    private:
        size_t n_sparse_points;
        std::map<SpeciesSet, TSpeciesGroup> species_groups;
        std::map<SpeciesSet, size_t> species_group_sizes;
    };

    template<typename TKernel>
    requires std::derived_from<TKernel, Kernel<TKernel::Dim, TKernel::Dependencies>>
    GenericKernelGroup<TKernel>::GenericKernelGroup(
        std::map<SpeciesSet, TSpeciesGroup> species_groups): species_groups(std::move(species_groups)) {

        n_sparse_points = 0;

        for (const auto& [species, hp_groups]: this->species_groups) {
            species_group_sizes[species] = 0;
            for (const auto& [hyper_params, sparse_points]: hp_groups) {
                species_group_sizes[species] += sparse_points.size();
            }
            n_sparse_points += species_group_sizes[species];
        }
    }

    template<typename TKernel>
    requires std::derived_from<TKernel, Kernel<TKernel::Dim, TKernel::Dependencies>>
    std::vector<AtomicQuantity> GenericKernelGroup<TKernel>::covariate(
        const size_t n_atoms, const std::map<SpeciesSet, std::vector<TDescriptor>> &descriptors) const {

        std::vector results(n_sparse_points, AtomicQuantity());

        size_t results_idx = 0;
        for (const auto& [species, hp_groups]: species_groups) {
            if (!descriptors.contains(species)) {
                results_idx += species_group_sizes[species];
                continue;
            }

            for (size_t tmp_idx = results_idx; tmp_idx < results_idx + species_group_sizes.at(species); tmp_idx++) {
                results[tmp_idx].forces.resize(n_atoms);
            }

            const auto& species_descriptors = descriptors.at(species);

            for (const auto& [parametrized_kernel, sparse_points]: hp_groups) {
                for (const auto& sparse_point: sparse_points) {
                    parametrized_kernel.covariance(results[results_idx++], sparse_point, species_descriptors);
                }
            }
        }

        return results;
    }

    template<typename TKernel>
    requires std::derived_from<TKernel, Kernel<TKernel::Dim, TKernel::Dependencies>>
    std::vector<Matrix> GenericKernelGroup<TKernel>::sparseToSparseCovariance() const {

        std::vector<Matrix> result;

        for (const auto& [species, hp_groups]: species_groups) {
            for (const auto& [parametrized_kernel, sparse_points]: hp_groups) {

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
        }

        return result;
    }
}

#endif