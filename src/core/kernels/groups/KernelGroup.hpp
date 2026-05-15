#ifndef JGAP_KERNELGROUP_HPP
#define JGAP_KERNELGROUP_HPP

#include "core/Matrix.hpp"
#include "core/atomic/Descriptor.hpp"

namespace jgap {

    template<size_t DescriptorDim, size_t DescriptorDependencies>
    class KernelGroup {
    public:
        using TDescriptor = Descriptor<DescriptorDim, DescriptorDependencies, CalculationType::WithGradients>;

        virtual ~KernelGroup() = default;

        virtual std::vector<AtomicQuantity> covariate(
            size_t n_atoms, const std::map<SpeciesSet, std::vector<TDescriptor>> &descriptors) const = 0;

        virtual std::vector<Matrix> sparseToSparseCovariance() const = 0;

        virtual size_t nSparsePoints() const = 0;
    };
}

#endif
