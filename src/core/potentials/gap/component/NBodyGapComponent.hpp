#ifndef JGAP_PERSPECIESKERNELGROUP_HPP
#define JGAP_PERSPECIESKERNELGROUP_HPP

#include "core/atomic/energy/AtomicQuantity.hpp"
#include "core/kernels/Kernel.hpp"
#include "core/transform/ClusterTransformation.hpp"
#include "GapComponent.hpp"
#include "core/Matrix.hpp"

namespace jgap {

    template<typename TTransformation, typename TKernel>
    requires CClusterTransformation<TTransformation> && CKernel<TKernel>
        && (TTransformation::Dim == TKernel::Dim)
    class NBodyGapComponent : public GapComponent {
    public:
        static constexpr size_t Dim = TKernel::Dim;
        static constexpr size_t ClusterSize = TTransformation::ClusterSize;
        static constexpr size_t Dependencies = Cluster<ClusterSize>::NSeparations;

        NBodyGapComponent(const SpeciesSet species, TTransformation transformer, TKernel kernel,
                          std::vector<Descriptor<Dim>> sparse_points,
                          const std::vector<Real> &optional_coeffs = {})
            : species(species), transformation(transformer), kernel(kernel), sparse_points(sparse_points) {

            if (!optional_coeffs.empty()) {
                setCoefficients(optional_coeffs);
            }
        }

        std::optional<AtomicQuantities> covariate(const NeighbourList &neighbour_list) const override {
            assert(neighbour_list.cutoff == getCutoff());

            auto clusters = neighbour_list.findAllClusters<ClusterSize>(species);

            if (clusters.empty()) {
                return std::nullopt;
            }

            AtomicQuantities result(nSparsePoints(), neighbour_list.nAtoms());

            for (const Cluster<ClusterSize>& cluster: clusters) {
                DescriptorAndDerivatives<Dim, ClusterSize> descriptor = transformation.evaluateAndDifferentiate(cluster);

                for (size_t sparse_idx = 0; sparse_idx < sparse_points.size(); sparse_idx++) {

                    const auto [K, gradK_wrt_q] = kernel.valueAndGradient(
                        sparse_points[sparse_idx].value, descriptor.value
                        );

                    result.energy(sparse_idx) += K;

                    for (size_t i = 0; i < Dependencies; i++) {
                        for (size_t j = i + 1; j < Dependencies; j++) {

                            const auto sep_idx = flattened_index(i, j);
                            const auto& separation = cluster.between(i, j);
                            const auto& derivatives_wrt_r_norms = descriptor.derivatives[sep_idx];

                            Real dK_drij_norm = 0.0;
                            for (size_t dim = 0; dim < Dim; dim++) {
                                dK_drij_norm += derivatives_wrt_r_norms[dim] * gradK_wrt_q[dim];
                            }

                            result.force(sparse_idx, cluster.atom_indexes[i]) -= separation.direction * dK_drij_norm;
                            result.force(sparse_idx, cluster.atom_indexes[j]) += separation.direction * dK_drij_norm;

                            result.virials(sparse_idx) += separation.virials * dK_drij_norm;
                        }
                    }
                }
            }

            return result;
        }

        Matrix sparseToSparseCovariance() const override {
            Matrix result(sparse_points.size(), sparse_points.size());

            for (size_t i = 0; i < sparse_points.size(); i++) {
                for (size_t j = i; j < sparse_points.size(); j++) {
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
            return transformation.getCutoffs();
        }

    private:
        SpeciesSet species;

        TTransformation transformation;
        TKernel kernel;

        std::vector<Descriptor<Dim>> sparse_points;
    };

    /*
    template<typename TKernel> requires CKernel<TKernel>
    class NBodyGapComponent : public SparseGroup<TKernel::Dim, TKernel::Dependencies> {
    public:
        constexpr static size_t DescriptorDim = TKernel::Dim;
        constexpr static size_t DescriptorDependencies = TKernel::Dependencies;

        using TDescriptor = Descriptor<DescriptorDim, DescriptorDependencies, CalculationType::WithGradients>;

        using TSparsePoint = std::array<Real, DescriptorDim>;
        using THyperParamGroup = std::pair<TKernel, std::vector<TSparsePoint>>;

        NBodyGapComponent(std::map<SpeciesSet, THyperParamGroup> species_groups);

        template<typename TTransformerGroup>
        requires CTransformerGroupOfDim<TTransformerGroup, TKernel::Dim, TKernel::Dependencies>
        static NBodyGapComponent constructFromTrainingData(
            const TTransformerGroup& transformation,
            const std::vector<NeighbourList>& training_data,
            const std::map<SpeciesSet, const Sparsifier<DescriptorDim>*>& sparsifiers = {},
            const Sparsifier<DescriptorDim>* default_sparsifier = nullptr,
            const std::map<SpeciesSet, TKernel>& kernels = {},
            const std::optional<TKernel>& default_kernel = std::nullopt,
            bool ignore_missing = false,
            bool ignore_two_body_symmetry = false);

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

    template<typename TKernel> requires CKernel<TKernel>
    NBodyGapComponent<TKernel>::NBodyGapComponent(
        std::map<SpeciesSet, THyperParamGroup> species_groups): species_groups(std::move(species_groups)) {

        n_sparse_points = 0;

        for (const auto& [species, hp_group]: this->species_groups) {
            species_group_sizes[species] = hp_group.second.size();
            n_sparse_points += species_group_sizes[species];
        }
    }

    template<typename TKernel> requires CKernel<TKernel>
    template<typename TTransformerGroup>
    requires CTransformerGroupOfDim<TTransformerGroup, TKernel::Dim, TKernel::Dependencies>
    NBodyGapComponent<TKernel> NBodyGapComponent<TKernel>::constructFromTrainingData(
        const TTransformerGroup& transformer_group,
        const std::vector<NeighbourList>& training_data,
        const std::map<SpeciesSet, const Sparsifier<DescriptorDim>*>& sparsifiers,
        const Sparsifier<DescriptorDim>* default_sparsifier,
        const std::map<SpeciesSet, TKernel>& kernels,
        const std::optional<TKernel>& default_kernel,
        const bool ignore_missing,
        const bool ignore_two_body_symmetry) {

        std::map<SpeciesSet, std::vector<Descriptor<DescriptorDim>>> all_descriptors;

        for (const auto& nl : training_data) {
            auto batch_descriptors = transformer_group.template transform<CalculationType::ValueOnly>(nl);

            for (const auto& [species_set, descriptors] : batch_descriptors) {
                if (species_set.isRedundantTwoBody() && !ignore_two_body_symmetry) {
                    continue;
                }

                for (const auto& desc : descriptors) {
                    Descriptor<DescriptorDim> flat_desc;
                    flat_desc.value = desc.value;
                    all_descriptors[species_set].push_back(flat_desc);
                }
            }
        }

        std::map<SpeciesSet, THyperParamGroup> species_groups;

        for (const auto& [species, descriptors] : all_descriptors) {

            const Sparsifier<DescriptorDim>* sparsifier_to_use = nullptr;
            if (sparsifiers.contains(species)) {
                sparsifier_to_use = sparsifiers.at(species);
            } else if (default_sparsifier != nullptr) {
                sparsifier_to_use = default_sparsifier;
            } else if (ignore_missing) {
                continue;
            } else {
                JGAP_LOG_AND_THROW("No sparsifier specified for species set '{}' and no default sparsifier provided.",
                    species.toString());
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
                JGAP_LOG_AND_THROW("No kernel specified for species set '{}' and no default kernel provided.",
                    species.toString());
            }

            species_groups[species] = { *kernel_to_use, std::move(sparse_points) };
        }

        return PerSpeciesKernelGroup(std::move(species_groups));
    }

    template<typename TKernel> requires CKernel<TKernel>
    std::vector<AtomicQuantity> NBodyGapComponent<TKernel>::covariate(
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

    template<typename TKernel> requires CKernel<TKernel>
    std::vector<Matrix> NBodyGapComponent<TKernel>::sparseToSparseCovariance() const {

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
    }*/
}

#endif