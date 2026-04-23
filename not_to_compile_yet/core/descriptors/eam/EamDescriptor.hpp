#ifndef EAMDESCRIPTOR_HPP
#define EAMDESCRIPTOR_HPP

#include <utility>
#include <queue>

#include "core/cutoff/CutoffFunction.hpp"
#include "core/descriptors/Descriptor.hpp"
#include "../../kernels/Kernel.hpp"
#include "data/descriptors/kernels/EamKernelIndex.hpp"
#include "pair_functions/EamPairFunction.hpp"
#include "io/parse/ParserRegistry.hpp"
#include "../../kernels/eam/EamSE.hpp"
#include "memory/MatrixBlock.hpp"

namespace jgap {

    class EamDescriptor : public Descriptor, Serializable {
    public:
        SETUP_PARSER_AND_SERIALIZATION(Descriptor, EamDescriptor, eam)

        EamDescriptor(std::vector<std::shared_ptr<EamKernel>> kernels,
                      std::shared_ptr<EamPairFunction> &default_pair_function,
                      std::map<ContributorReceiverSpecies,
                      std::shared_ptr<EamPairFunction>> pairFunctions);

        CutoffRanges getCutoff() override;;

        std::vector<std::shared_ptr<IKernel>> getKernels() override;
        void setupSparseKernels(const std::vector<AtomicStructure> &fromData) override;

        std::vector<Covariance> covariate(const AtomicStructure &atomicStructure) override;
        std::vector<std::shared_ptr<MatrixBlock>> selfCovariate() override;

        Predictions predict(const AtomicStructure &atomicStructure) override;

        void tabulate(TabulationData &table) override;

    private:
        double max_cutoff_;
        std::vector<std::shared_ptr<EamKernel>> kernels_;
        std::map<Species, std::vector<size_t>> kernel_indices_per_species_;

        // A queue to emphasize the one-time use
        std::queue<DataNode> _kernelSetups;

        std::shared_ptr<EamPairFunction> default_pair_function_;
        std::map<ContributorReceiverSpecies, std::shared_ptr<EamPairFunction>> pair_functions_;

        EamKernelIndex doIndex(const AtomicStructure &structure) const;

        void mapKernelIds();
    };
}



#endif
