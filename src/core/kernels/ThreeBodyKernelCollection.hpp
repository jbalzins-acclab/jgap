#ifndef JGAP_THREEBODYKERNELCOLLECTION_HPP
#define JGAP_THREEBODYKERNELCOLLECTION_HPP
#include "KernelCollection.hpp"
#include "core/cutoff/CutoffFunction.hpp"
#include "core/descriptors/3b/ThreeBodyDescriptorFinder.hpp"

namespace jgap {

    class ThreeBodyKernelCollection : public KernelCollection/*, Tabulatable, Serializable*/ {
    public:

        ThreeBodyKernelCollection(
            std::shared_ptr<ThreeBodyDescriptorFinder> descriptor_transformer,
            std::map<EncodedSpeciesSets, std::vector<std::shared_ptr<RealKernel<3, 3>>>> kernels_per_species_triplet
            )
            : descriptor_finder_(std::move(descriptor_transformer)),
              kernels_per_species_triplet_(std::move(kernels_per_species_triplet)) {
            //optional_angular_cutoff_(std::move(optional_angular_cutoff)) {

        }

        std::vector<Predictions> covariate(const AtomicStructure &atomic_structure) override {

            auto qs_per_species = descriptor_finder_->find(atomic_structure);

            std::vector<Predictions> result;
            for (const auto& [species_triplet, descriptors]: qs_per_species) {
                if (!kernels_per_species_triplet_.contains(species_triplet)) {
                    continue;
                }
                for (const auto& kernel: kernels_per_species_triplet_[species_triplet]) {
                    result.push_back(kernel->covariance(atomic_structure, descriptors));
                }
            }
            return result;
        }

        std::vector<std::shared_ptr<MatrixBlock>> selfCovariate() override {
            std::vector<std::shared_ptr<MatrixBlock>> result;
            for (const auto& kernels_per_species: kernels_per_species_triplet_ | std::views::values) {
                size_t n = kernels_per_species.size();

                auto block = std::make_shared<MatrixBlock>(n, n);

                for (size_t i = 0; i < n; i++) {
                    for (size_t j = i; j < n; j++) {
                        (*block)(i, j) = kernels_per_species[i]->crossCovariance(kernels_per_species[j]);
                    }
                }

                result.push_back(block);
            }
            return result;
        }

        std::vector<std::shared_ptr<IKernel>> getKernels() override {
            std::vector<std::shared_ptr<IKernel>> result;

            for (const auto& kernels_per_species: kernels_per_species_triplet_ | std::views::values) {
                for (const auto& kernel: kernels_per_species) {
                    result.push_back(std::static_pointer_cast<IKernel>(kernel));
                }
            }

            return result;
        }

    private:
        std::shared_ptr<ThreeBodyDescriptorFinder> descriptor_finder_;
        std::map<EncodedSpeciesSets, std::vector<std::shared_ptr<RealKernel<3, 3>>>> kernels_per_species_triplet_;
    };
}

#endif
