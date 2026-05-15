#ifndef JGAP_BLINDKERNELGROUP_HPP
#define JGAP_BLINDKERNELGROUP_HPP

#include "core/atomic/energy/AtomicQuantity.hpp"
#include "core/kernels/Kernel.hpp"
#include "core/transform/Transformer.hpp"
#include "core/sparsification/Sparsifier.hpp"
#include <concepts>

#include "KernelGroup.hpp"
#include "core/Matrix.hpp"

namespace jgap {

    template<typename TKernel>
    requires std::derived_from<TKernel, Kernel<TKernel::Dim, TKernel::Dependencies>>
    class BlindKernelGroup : public KernelGroup<TKernel::Dim, TKernel::Dependencies> {
    public:
        constexpr static size_t DescriptorDim = TKernel::Dim;
        constexpr static size_t DescriptorDependencies = TKernel::Dependencies;

        using TDescriptor = Descriptor<DescriptorDim, DescriptorDependencies, CalculationType::WithGradients>;

        using TSparsePoint = std::array<Real, DescriptorDim>;

        BlindKernelGroup(TKernel kernel, std::vector<TSparsePoint> sparse_points);

        static BlindKernelGroup constructFromTrainingData(
            TKernel kernel,
            const std::shared_ptr<Transformer<DescriptorDim, DescriptorDependencies>>& transformer,
            const std::shared_ptr<Sparsifier<DescriptorDim>>& sparsifier,
            const std::vector<NeighbourList>& training_data);

        std::vector<AtomicQuantity> covariate(size_t n_atoms,
            const std::map<SpeciesSet, std::vector<TDescriptor>> &descriptors) const override;
        std::vector<Matrix> sparseToSparseCovariance() const override;

        size_t nSparsePoints() const override {
            return n_sparse_points;
        }

    private:
        size_t n_sparse_points;
        TKernel kernel;
        std::vector<TSparsePoint> sparse_points;
    };

    template<typename TKernel>
    requires std::derived_from<TKernel, Kernel<TKernel::Dim, TKernel::Dependencies>>
    BlindKernelGroup<TKernel>::BlindKernelGroup(
        TKernel kernel, std::vector<TSparsePoint> sparse_points):
        kernel(std::move(kernel)), sparse_points(std::move(sparse_points)) {

        n_sparse_points = this->sparse_points.size();
    }

    template<typename TKernel>
    requires std::derived_from<TKernel, Kernel<TKernel::Dim, TKernel::Dependencies>>
    BlindKernelGroup<TKernel> BlindKernelGroup<TKernel>::constructFromTrainingData(
        TKernel kernel,
        const std::shared_ptr<Transformer<DescriptorDim, DescriptorDependencies>>& transformer,
        const std::shared_ptr<Sparsifier<DescriptorDim>>& sparsifier,
        const std::vector<NeighbourList>& training_data) {

        using TValueOnly = Transformer<DescriptorDim, DescriptorDependencies>::TValueOnly;
        std::vector<Descriptor<DescriptorDim>> all_descriptors;

        for (const auto& nl : training_data) {
            auto batch_descriptors = transformer->transform(nl);
            for (const auto& [species, descriptors] : batch_descriptors) {
                for (const auto& desc : descriptors) {
                    Descriptor<DescriptorDim> flat_desc;
                    flat_desc.value = desc.value;
                    all_descriptors.push_back(flat_desc);
                }
            }
        }

        auto sparse_descriptors = sparsifier->selectSparsePoints(all_descriptors);
        std::vector<TSparsePoint> sparse_points;
        sparse_points.reserve(sparse_descriptors.size());
        for (const auto& desc : sparse_descriptors) {
            sparse_points.push_back(desc.value);
        }

        return BlindKernelGroup(std::move(kernel), std::move(sparse_points));
    }

    template<typename TKernel>
    requires std::derived_from<TKernel, Kernel<TKernel::Dim, TKernel::Dependencies>>
    std::vector<AtomicQuantity> BlindKernelGroup<TKernel>::covariate(
        const size_t n_atoms, const std::map<SpeciesSet, std::vector<TDescriptor>> &descriptors) const {

        std::vector results(n_sparse_points, AtomicQuantity(n_atoms));

        for (size_t tmp_index = 0; tmp_index < n_sparse_points; tmp_index++) {
            results[tmp_index].forces.resize(n_atoms);
        }

        size_t results_idx = 0;
        for (const auto& sparse_point: sparse_points) {
            for (const auto& [species, species_descriptors]: descriptors) {
                kernel.covariance(results[results_idx], sparse_point, species_descriptors);
            }
            results_idx++;
        }

        return results;
    }

    template<typename TKernel>
    requires std::derived_from<TKernel, Kernel<TKernel::Dim, TKernel::Dependencies>>
    std::vector<Matrix> BlindKernelGroup<TKernel>::sparseToSparseCovariance() const {

        std::vector<Matrix> result;

        result.emplace_back(sparse_points.size(), sparse_points.size());
        auto& current_matrix = result.back();
        size_t m = sparse_points.size();

        for (size_t i = 0; i < m; i++) {
            for (size_t j = i; j < m; j++) {
                current_matrix(i, j) = kernel.value(sparse_points[i], sparse_points[j]);
                current_matrix(j, i) = current_matrix(i, j);
            }
        }

        return result;
    }
}

#endif // JGAP_BLINDKERNELGROUP_HPP
