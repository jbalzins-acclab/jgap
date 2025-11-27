#ifndef EAMDESCRIPTOR_HPP
#define EAMDESCRIPTOR_HPP

#include <nlohmann/json.hpp>
#include <utility>
#include <queue>

#include "core/cutoff/CutoffFunction.hpp"
#include "core/descriptors/Descriptor.hpp"
#include "../Kernel.hpp"
#include "data/descriptors/kernels/EamKernelIndex.hpp"
#include "pair_functions/EamPairFunction.hpp"
#include "io/parse/ParserRegistry.hpp"
#include "EamSE.hpp"
#include "memory/MatrixBlock.hpp"
#include "../sparsification/Sparsifier.hpp"

namespace jgap {

    class EamDescriptor : public Descriptor {
    public:
        static constexpr string TYPE = "eam";
        EamDescriptor(vector<shared_ptr<EamKernel>> kernels,
                      shared_ptr<EamPairFunction> &defaultPairFunction,
                      map<ContributorReceiverSpecies, shared_ptr<EamPairFunction>> pairFunctions);

        EamDescriptor(const nlohmann::json &params);
        nlohmann::json serialize() override;
        string getType() override { return TYPE; }

        CutoffRanges getCutoff() override;;

        vector<shared_ptr<IKernel>> getKernels() override;
        void setupSparseKernels(const vector<AtomicStructure> &fromData) override;

        vector<Covariance> covariate(const AtomicStructure &atomicStructure) override;
        vector<shared_ptr<MatrixBlock>> selfCovariate() override;

        Predictions predict(const AtomicStructure &atomicStructure) override;

        void tabulate(TabulationData &table) override;

    private:
        double _maxCutoff;
        vector<shared_ptr<EamKernel>> _kernels;
        map<Species, vector<size_t>> _kernelIndicesPerSpecies;

        // A queue to emphasize the one-time use
        queue<nlohmann::json> _kernelSetups;

        shared_ptr<EamPairFunction> _defaultPairFunction;
        map<ContributorReceiverSpecies, shared_ptr<EamPairFunction>> _pairFunctions;

        EamKernelIndex doIndex(const AtomicStructure &structure) const;

        void mapKernelIds();
    };

    REGISTER_PARSER(Descriptor, EamDescriptor)
}



#endif
