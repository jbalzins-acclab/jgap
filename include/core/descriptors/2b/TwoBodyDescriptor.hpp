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
#include <nlohmann/json.hpp>

namespace jgap {

    class TwoBodyDescriptor : public Descriptor {
    public:
        static constexpr string TYPE = "2b";
        TwoBodyDescriptor(shared_ptr<CutoffFunction> cutoffFunction, vector<shared_ptr<TwoBodyKernel>> kernels);
        TwoBodyDescriptor(const nlohmann::json& params);
        nlohmann::json serialize() override;
        string getType() override { return TYPE; };

        CutoffRanges getCutoff() override { return CutoffRanges{.twoBody = _cutoffFunction->getCutoff()}; }

        vector<shared_ptr<IKernel>> getKernels() override;
        void setupKernels(const vector<AtomicStructure> &fromData) override;

        vector<Covariance> covariate(const AtomicStructure &atomicStructure) override;
        vector<shared_ptr<MatrixBlock>> selfCovariate() override;

        Predictions predict(const AtomicStructure &atomicStructure) override;

        void tabulate(TabulationData &table) override;

    protected:
        shared_ptr<CutoffFunction> _cutoffFunction;
        vector<shared_ptr<TwoBodyKernel>> _kernels;
        map<SpeciesPair, vector<size_t>> _kernelIdsPerSpeciesPair;

        // A queue to emphasize the one-time use
        queue<nlohmann::json> _kernelSetups;

        map<SpeciesPair, TwoBodyIndex> doIndex(const AtomicStructure &atomicStructure) const;
        bool checkSpecies(SpeciesPair pairInData, nlohmann::json filters);

        void mapKernelIds();
    };

    REGISTER_PARSER(Descriptor, TwoBodyDescriptor)
}

#endif //TWOBODYDESCRIPTOR_HPP
