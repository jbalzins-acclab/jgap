#ifndef THREEBODYDESCRIPTOR_HPP
#define THREEBODYDESCRIPTOR_HPP

#include <nlohmann/json.hpp>
#include <utility>

#include "core/descriptors/Descriptor.hpp"
#include "core/kernels/ThreeBodySE.hpp"

#include "data/params/ThreeBodyDescriptorParams.hpp"
#include "memory/MatrixBlock.hpp"

namespace jgap {

    class ThreeBodyDescriptor : public Descriptor {
    public:
        explicit ThreeBodyDescriptor(ThreeBodyDescriptorParams params);
        ~ThreeBodyDescriptor() override = default;

        double getCutoff() override { return _params.cutoff; }

        size_t nSparsePoints() override;

        nlohmann::json serialize() override;

        void setSparsePoints(const vector<AtomicStructure>& fromData) override;

        void setSparsePoints(map<SpeciesTriplet, vector<Vector3>> sparsePointsPerSpeciesPair) {
            _sparsePointsPerSpeciesTriplet = std::move(sparsePointsPerSpeciesPair);
        }

        vector<Covariance> covariate(const AtomicStructure &atomicStructure) override;
        vector<pair<size_t, shared_ptr<MatrixBlock>>> selfCovariate() override;

    private:
        string _name;
        ThreeBodyDescriptorParams _params;
        shared_ptr<Kernel<Vector3, ThreeBodyKernelIndex>> _kernel;

        map<SpeciesTriplet, vector<Vector3>> _sparsePointsPerSpeciesTriplet;

        void sparsifyTrueUniform(const map<SpeciesTriplet, vector<Vector3>> &all3b);
        void sparsifyQuipUniform(const map<SpeciesTriplet, vector<Vector3>>& all3b);

        [[nodiscard]]
        map<SpeciesTriplet, ThreeBodyKernelIndex> doIndex(const AtomicStructure &atomicStructure) const;
    };
}

#endif //THREEBODYDESCRIPTOR_HPP
