#ifndef THREEBODYDESCRIPTOR_HPP
#define THREEBODYDESCRIPTOR_HPP

#include <nlohmann/json.hpp>
#include <utility>

#include "core/descriptors/Descriptor.hpp"
#include "core/descriptors/kernels/ThreeBodySE.hpp"

#include "data/params/ThreeBodyDescriptorParams.hpp"
#include "io/parse/ParserRegistry.hpp"
#include "memory/MatrixBlock.hpp"
#include "sparsification/3b/PerSpecies3bSparsifier.hpp"

namespace jgap {

    class ThreeBodyDescriptor : public Descriptor {
    public:
        ~ThreeBodyDescriptor() override = default;

        explicit ThreeBodyDescriptor(const nlohmann::json& params);
        nlohmann::json serialize() override;
        string getType() override { return "3b"; }

        double getCutoff() override { return _cutoff; }

        size_t nSparsePoints() override;
        void setSparsePoints(const vector<AtomicStructure>& fromData) override;
        void setSparsePoints(map<SpeciesTriplet, vector<Vector3>> sparsePointsPerSpeciesPair) {
            _sparsePointsPerSpeciesTriplet = std::move(sparsePointsPerSpeciesPair);
        }

        vector<Covariance> covariate(const AtomicStructure &atomicStructure) override;
        vector<pair<size_t, shared_ptr<MatrixBlock>>> selfCovariate() override;

    private:
        double _cutoff;
        shared_ptr<Kernel<Vector3, ThreeBodyKernelIndex>> _kernel;
        shared_ptr<PerSpecies3bSparsifier> _sparsifier;

        map<SpeciesTriplet, vector<Vector3>> _sparsePointsPerSpeciesTriplet;

        [[nodiscard]]
        map<SpeciesTriplet, ThreeBodyKernelIndex> doIndex(const AtomicStructure &atomicStructure) const;
    };

    REGISTER_PARSER("3b", Descriptor, ThreeBodyDescriptor)
}

#endif //THREEBODYDESCRIPTOR_HPP
