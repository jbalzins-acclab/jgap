#ifndef JGAP_KERNELCOLLECTION_HPP
#define JGAP_KERNELCOLLECTION_HPP
#include <vector>

#include "Kernel.hpp"
#include "data/CutoffRanges.hpp"
#include "data/atomic/PredictionData.hpp"
#include "memory/MatrixBlock.hpp"

namespace jgap {


    class KernelCollection {
    public:
        virtual ~KernelCollection() = default;

        virtual std::vector<Predictions> covariate(const AtomicStructure &atomic_structure) = 0;
        virtual std::vector<std::shared_ptr<MatrixBlock>> selfCovariate() = 0;

        virtual CutoffRanges getCutoff();
        // Must be in same order as @covariate and @selfCovariate
        virtual std::vector<std::shared_ptr<IKernel>> getKernels() = 0;

        virtual Predictions predict(const AtomicStructure &atomic_structure) {

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
    };
}

#endif