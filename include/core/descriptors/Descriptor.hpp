#ifndef DESCRIPTOR_HPP
#define DESCRIPTOR_HPP

#include "data/Vector3.hpp"
#include "memory/MatrixBlock.hpp"

#include <memory>
#include <vector>
#include <nlohmann/json_fwd.hpp>

#include "Kernel.hpp"
#include "data/CutoffRanges.hpp"
#include "data/TabulationData.hpp"

using namespace std;

namespace jgap {

    class Descriptor {
    public:
        virtual ~Descriptor() = default;

        virtual nlohmann::json serialize() = 0;
        virtual string getType() = 0;
        virtual CutoffRanges getCutoff() = 0;
        // Must be in same order as @covariate and @selfCovariate
        virtual vector<shared_ptr<IKernel>> getKernels() = 0;

        // Sparsification strategy to constructor
        virtual void setupKernels(const vector<AtomicStructure> &fromData) = 0;

        virtual vector<Covariance> covariate(const AtomicStructure &atomicStructure) = 0;
        virtual vector<shared_ptr<MatrixBlock>> selfCovariate() = 0;

        virtual Predictions predict(const AtomicStructure &atomicStructure) = 0;

        virtual void tabulate(TabulationData &table) = 0;
    };
}

#endif
