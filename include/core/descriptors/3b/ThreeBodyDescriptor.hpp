#ifndef THREEBODYDESCRIPTOR_HPP
#define THREEBODYDESCRIPTOR_HPP

#include <nlohmann/json.hpp>
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
        static constexpr string TYPE = "3b";

        ThreeBodyDescriptor(shared_ptr<CutoffFunction>& cutoffFunction, vector<shared_ptr<ThreeBodyKernel>>& kernels);
        ThreeBodyDescriptor(const nlohmann::json& params);

        nlohmann::json serialize() override;
        string getType() override { return TYPE; }

        CutoffRanges getCutoff() override { return CutoffRanges{.threeBody = _cutoffFunction->getCutoff()}; };

        void setupKernels(const vector<AtomicStructure>& fromData) override;
        vector<shared_ptr<IKernel>> getKernels() override;

        Predictions predict(const AtomicStructure &atomicStructure) override;

        vector<Covariance> covariate(const AtomicStructure &atomicStructure) override;
        vector<shared_ptr<MatrixBlock>> selfCovariate() override;

        void tabulate(TabulationData &table) override;

        double invariantTripletToCutoff(const Vector3 &t) const;
        static Vector3 toInvariantTriplet(double r01, double r02, double r12);
        static array<Vector3, 3> invariantTripletGradients(double r01, double r02);

    protected:
        double _cutoff;
        shared_ptr<CutoffFunction> _cutoffFunction;
        vector<shared_ptr<ThreeBodyKernel>> _kernels;
        map<SpeciesTriplet, vector<size_t>> _kernelIdsPerSpeciesTriplet;

        shared_ptr<Sparsifier> _sparsifier;

        // A queue to emphasize the one-time use
        queue<nlohmann::json> _kernelSetups;

        map<SpeciesTriplet, ThreeBodyIndex> doIndex(const AtomicStructure &atomicStructure) const;
        bool checkSpecies(const SpeciesTriplet& tripletInData, nlohmann::json filters);

        void mapKernelIds();
    };

    REGISTER_PARSER(Descriptor, ThreeBodyDescriptor)
}

#endif
