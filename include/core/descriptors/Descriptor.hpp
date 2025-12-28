#ifndef DESCRIPTOR_HPP
#define DESCRIPTOR_HPP

#include "data/Vector3.hpp"
#include "memory/MatrixBlock.hpp"

#include <memory>
#include <vector>
#include <string>
#include "data/DataNode.hpp"

#include "Kernel.hpp"
#include "data/CutoffRanges.hpp"
#include "../../data/tabulation/TabulationData.hpp"

namespace jgap {

    class Descriptor {
    public:
        virtual ~Descriptor() = default;

        virtual CutoffRanges getCutoff() = 0;
        // Must be in same order as @covariate and @selfCovariate
        virtual std::vector<std::shared_ptr<IKernel>> getKernels() = 0;

        // Sparsification strategy to constructor
        virtual void setupSparseKernels(const std::vector<AtomicStructure> &from_data) = 0;

        virtual std::vector<Covariance> covariate(const AtomicStructure &atomic_structure) = 0;
        virtual std::vector<std::shared_ptr<MatrixBlock>> selfCovariate() = 0;

        virtual Predictions predict(const AtomicStructure &atomic_structure) = 0;

        virtual void tabulate(TabulationData &table) = 0;
    };
}

#endif
