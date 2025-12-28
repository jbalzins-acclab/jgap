#ifndef TWOBODYDESCRIPTOR_HPP
#define TWOBODYDESCRIPTOR_HPP

#include "core/descriptors/Descriptor.hpp"
#include "../Kernel.hpp"
#include "core/cutoff/CutoffFunction.hpp"
#include "data/descriptors/kernels/TwoBodyIndex.hpp"
#include "io/parse/ParserRegistry.hpp"
#include "TwoBodySE.hpp"
#include "memory/MatrixBlock.hpp"
#include "../sparsification/Sparsifier.hpp"

#include <ranges>
#include <utility>
#include <set>
#include <queue>

namespace jgap {

    class TwoBodyDescriptor : public Descriptor {
    public:
        static constexpr std::string TYPE = "2b";
        static std::shared_ptr<TwoBodyDescriptor> fromDataNode(const DataNode& params);

        TwoBodyDescriptor(std::shared_ptr<CutoffFunction> cutoffFunction, std::vector<std::shared_ptr<TwoBodyKernel>> kernels);

        DataNode serialize() override;

        std::string getType() override { return TYPE; };

        CutoffRanges getCutoff() override { return CutoffRanges{.twoBody = _cutoffFunction->getCutoff()}; }

        std::vector<std::shared_ptr<IKernel>> getKernels() override;
        void setupSparseKernels(const std::vector<AtomicStructure> &fromData) override;

        std::vector<Covariance> covariate(const AtomicStructure &atomicStructure) override;
        std::vector<std::shared_ptr<MatrixBlock>> selfCovariate() override;

        Predictions predict(const AtomicStructure &atomicStructure) override;

        void tabulate(TabulationData &table) override;

    protected:
        std::shared_ptr<CutoffFunction> _cutoffFunction;
        std::vector<std::shared_ptr<TwoBodyKernel>> _kernels;
        std::map<SpeciesPair, std::vector<size_t>> _kernelIdsPerSpeciesPair;

        // A queue to emphasize the one-time use
        std::queue<DataNode> _kernelSetups;

        std::map<SpeciesPair, TwoBodyIndex> doIndex(const AtomicStructure &atomicStructure) const;
        bool checkSpecies(const SpeciesPair& pairInData, DataNode filters);

        void mapKernelIds();
    };

    SETUP_PARSER(Descriptor, TwoBodyDescriptor)
}

#endif
