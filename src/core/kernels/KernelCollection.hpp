#ifndef JGAP_KERNELCOLLECTION_HPP
#define JGAP_KERNELCOLLECTION_HPP

#include <ranges>
#include <vector>
#include <concepts>
#include <type_traits>

#include "Kernel.hpp"
#include "core/descriptors/SeparationsFinder.hpp"
#include "data/CutoffRanges.hpp"
#include "data/atomic/Predictions.hpp"
#include "../Matrix.hpp"

namespace jgap {

    class IKernelCollection {
    public:
        virtual ~IKernelCollection() = default;
        virtual std::vector<std::shared_ptr<IKernel>> getKernels() = 0;
        virtual std::vector<Predictions> covariate(const Box &atomic_structure) = 0;
        virtual std::vector<std::shared_ptr<Matrix>> selfCovariate() = 0;
        virtual CutoffRanges getCutoff() = 0;
        virtual Predictions predict(const Box &atomic_structure) = 0;
    };

    template <size_t N_DIMENSIONS, typename T>
    concept DerivedFromFinder = std::derived_from<T, SeparationsFinder<N_DIMENSIONS>>;

    /**
     * Collection of kernels that share the same cutoff logic, and dimension.
     * @tparam N_DIMENSIONS
     * @tparam TDescriptorFinder
     */
    template<size_t N_DIMENSIONS, typename TDescriptorFinder>
    requires DerivedFromFinder<N_DIMENSIONS, TDescriptorFinder>
    class KernelCollection : public IKernelCollection {
    public:

        KernelCollection(
            std::shared_ptr<TDescriptorFinder> descriptor_transformer,
            std::map<EncodedSpeciesSets, std::vector<std::shared_ptr<RealKernel<N_DIMENSIONS>>>> kernels_per_species_set
            )
            : descriptor_finder_(std::move(descriptor_transformer)),
              kernels_per_species_set_(std::move(kernels_per_species_set)) {
            //optional_angular_cutoff_(std::move(optional_angular_cutoff)) {
            n_kernels_ = 0;
            for (const auto& kernels: kernels_per_species_set_ | std::views::values) {
                n_kernels_ += kernels.size();
            }
        }

        std::vector<Predictions> covariate(const Box &atomic_structure) override {

            auto descriptors_per_species = descriptor_finder_->find(atomic_structure);

            size_t index = 0;
            std::vector<Predictions> result(n_kernels_);
            for (const auto& [species_triplet, kernels]: kernels_per_species_set_) {
                if (!descriptors_per_species.contains(species_triplet)) {
                    index += kernels.size();
                    continue;
                }

                //JGAP_LOG_WARN("bbb");
                const auto& descriptors = descriptors_per_species[species_triplet];

                for (const auto& kernel: kernels) {
                    result[index++] = kernel->covariance(atomic_structure, descriptors);
                }
            }
            //JGAP_LOG_WARN("aaa{}",result.size());
            return result;
        }

        std::vector<std::shared_ptr<Matrix> > selfCovariate() override {

            std::vector<std::shared_ptr<Matrix>> result;
            for (const auto& kernels_per_species: kernels_per_species_set_ | std::views::values) {
                size_t n = kernels_per_species.size();

                auto block = std::make_shared<Matrix>(n, n);

                for (size_t i = 0; i < n; i++) {
                    for (size_t j = i; j < n; j++) {
                        (*block)(i, j) = kernels_per_species[i]->crossCovariance(kernels_per_species[j]);
                    }
                }

                result.push_back(block);
            }

            return result;
        }

        CutoffRanges getCutoff() override {
            return descriptor_finder_->getCutoff();
        }

        // Must be in same order as @covariate and @selfCovariate
        std::vector<std::shared_ptr<IKernel>> getKernels() override {
            std::vector<std::shared_ptr<IKernel>> result;

            for (const auto& kernels_per_species: kernels_per_species_set_ | std::views::values) {
                for (const auto& kernel: kernels_per_species) {
                    result.push_back(std::static_pointer_cast<IKernel>(kernel));
                }
            }

            return result;
        }

        Predictions predict(const Box &atomic_structure) override {

            const auto kernels = getKernels();
            const auto covariances = covariate(atomic_structure);

            assert(kernels.size() == covariances.size() && "Kernel - covariance size mismatch");
            const size_t n_kernels = kernels.size();

            Predictions result{};
            for (size_t i = 0; i < n_kernels; i++) {
                assert(kernels[i]->coefficient.has_value() && "Kernel coefficient not set");
                result += covariances[i] * kernels[i]->coefficient.value();
            }

            return result;
        }

    protected:
        std::shared_ptr<TDescriptorFinder> descriptor_finder_;
        std::map<EncodedSpeciesSets, std::vector<std::shared_ptr<RealKernel<N_DIMENSIONS>>>> kernels_per_species_set_;
        size_t n_kernels_;
    };
}

#endif