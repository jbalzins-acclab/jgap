#ifndef THREEBODYDESCRIPTOR_HPP
#define THREEBODYDESCRIPTOR_HPP

#include <nlohmann/json.hpp>
#include <utility>

#include "core/descriptors/Descriptor.hpp"
#include "core/descriptors/kernels/ThreeBodySE.hpp"
#include "data/descriptors/kernels/ThreeBodyKernelIndex.hpp"

#include "io/parse/ParserRegistry.hpp"
#include "memory/MatrixBlock.hpp"
#include "sparsification/Sparsifier.hpp"

namespace jgap {

    class ThreeBodyDescriptor : public Descriptor {
    public:
        explicit ThreeBodyDescriptor(const nlohmann::json& params);
        nlohmann::json serialize() override;
        string getType() override { return "3b"; }

        double getCutoff() override { return _cutoffFunction->getCutoff(); }

        size_t nSparsePoints() override;
        void setSparsePoints(const vector<AtomicStructure>& fromData) override;
        void setSparsePoints(const map<SpeciesTriplet, vector<Vector3>>& sparsePointsPerSpeciesTriplet);

        vector<Covariance> covariate(const AtomicStructure &atomicStructure) override;
        vector<pair<size_t, shared_ptr<MatrixBlock>>> selfCovariate() override;

        TabulationData tabulate(const TabulationParams &params) override;

        static Vector3 toInvariantTriplet(double r01, double r02, double r12);
        static array<Vector3, 3> invariantTripletGradients(double r01, double r02);
        double invariantTripletToCutoff(const Vector3 &t) const;

    private:
        double _cutoff;
        shared_ptr<CutoffFunction> _cutoffFunction;
        shared_ptr<ThreeBodyKernel> _kernel;
        shared_ptr<Sparsifier> _sparsifier;

        map<SpeciesTriplet, vector<ThreeBodyDescriptorData>> _sparsePointsPerSpeciesTriplet;

        [[nodiscard]]
        map<SpeciesTriplet, ThreeBodyKernelIndex> doIndex(const AtomicStructure &atomicStructure) const;
    };

    REGISTER_PARSER("3b", Descriptor, ThreeBodyDescriptor)
}

#endif //THREEBODYDESCRIPTOR_HPP
