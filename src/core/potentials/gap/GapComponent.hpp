#ifndef JGAP_GAPCOMPONENT_HPP
#define JGAP_GAPCOMPONENT_HPP

#include <vector>

#include "core/atomic/Atoms.hpp"
#include "core/atomic/energy/AtomicQuantity.hpp"
#include "core/atomic/neighbours/NeighbourList.hpp"
#include "../../kernels/groups/KernelGroup.hpp"
#include "core/transform/Transformer.hpp"
#include "io/Serializable.hpp"

namespace jgap {

    class IGapComponent {
    public:
        using Ptr = std::shared_ptr<IGapComponent>;

        virtual ~IGapComponent() = default;
        virtual std::vector<AtomicQuantity> covariate(const NeighbourList& neighbour_list) const = 0;
        virtual std::vector<Matrix> sparseToSparseCovariance() const = 0;
        virtual size_t nSparsePoints() const = 0;
        virtual Cutoffs getCutoffs() const = 0;
    };

    template<size_t DescriptorDim, size_t DescriptorDependencies>
    class GapComponent : public IGapComponent {
    public:
        //static constexpr size_t DescriptorDim = 4;
        //static constexpr size_t DescriptorDependencies = 3;
        using TTransformer = Transformer<DescriptorDim, DescriptorDependencies>;
        using TKernelGroup = KernelGroup<DescriptorDim, DescriptorDependencies>;

        GapComponent(std::shared_ptr<TTransformer> transformer, std::shared_ptr<TKernelGroup> kernel_group)
            : transformer(std::move(transformer)), kernel_group(std::move(kernel_group)) {
        }

        size_t nSparsePoints() const override {
            return kernel_group->nSparsePoints();
        }

        std::vector<Matrix> sparseToSparseCovariance() const override {
            return kernel_group->sparseToSparseCovariance();
        }

        Cutoffs getCutoffs() const override {
            return transformer->getCutoffs();
        }

        std::vector<AtomicQuantity> covariate(const NeighbourList& neighbour_list) const override {
            auto descriptors_per_species = transformer->transformWithGradients(neighbour_list);
            return kernel_group->covariate(neighbour_list.nAtoms(), descriptors_per_species);
        }

    private:
        std::shared_ptr<TTransformer> transformer;
        std::shared_ptr<TKernelGroup> kernel_group;
    };
}

#endif
