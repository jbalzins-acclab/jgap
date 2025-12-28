#ifndef THREEBODYDESCRIPTOR_HPP
#define THREEBODYDESCRIPTOR_HPP

#include <queue>
#include <utility>

#include "core/descriptors/Descriptor.hpp"
#include "ThreeBodySE.hpp"
#include "data/descriptors/kernels/ThreeBodyIndex.hpp"

#include "io/parse/ParserRegistry.hpp"
#include "memory/MatrixBlock.hpp"
#include "../sparsification/Sparsifier.hpp"

namespace jgap {

    class ThreeBodyDescriptor : public Descriptor {
    public:
        static constexpr std::string TYPE = "3b";
        static std::shared_ptr<ThreeBodyDescriptor> fromDataNode(const DataNode& params);

        ThreeBodyDescriptor(std::shared_ptr<CutoffFunction> cutoffFunction,
                            std::vector<std::shared_ptr<ThreeBodyKernel>>& kernels);

        DataNode serialize() override;
        std::string getType() override { return TYPE; }

        CutoffRanges getCutoff() override { return CutoffRanges{.threeBody = _cutoffFunction->getCutoff()}; };

        void setupSparseKernels(const std::vector<AtomicStructure>& fromData) override;
        std::vector<std::shared_ptr<IKernel>> getKernels() override;

        Predictions predict(const AtomicStructure &atomicStructure) override;

        std::vector<Covariance> covariate(const AtomicStructure &atomicStructure) override;
        std::vector<std::shared_ptr<MatrixBlock>> selfCovariate() override;

        void tabulate(TabulationData &table) override;

        double invariantTripletToCutoff(const Vector3 &t) const;
        static Vector3 toInvariantTriplet(double r01, double r02, double r12);
        static std::array<Vector3, 3> invariantTripletGradients(double r01, double r02);

    protected:
        double _cutoff;
        std::shared_ptr<CutoffFunction> _cutoffFunction;
        std::vector<std::shared_ptr<ThreeBodyKernel>> _kernels;
        std::map<SpeciesTriplet, std::vector<size_t>> _kernelIdsPerSpeciesTriplet;

        std::shared_ptr<Sparsifier> _sparsifier;

        // A queue to emphasize the one-time use
        std::queue<DataNode> _kernelSetups;

        std::map<SpeciesTriplet, ThreeBodyIndex> doIndex(const AtomicStructure &atomicStructure) const;
        bool checkSpecies(const SpeciesTriplet& tripletInData, DataNode filters);

        void mapKernelIds();
    };

    SETUP_PARSER(Descriptor, ThreeBodyDescriptor)
}

#endif
