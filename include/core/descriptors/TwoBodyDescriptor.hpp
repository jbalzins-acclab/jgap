#ifndef TWOBODYDESCRIPTOR_HPP
#define TWOBODYDESCRIPTOR_HPP

#include <ranges>
#include <utility>
#include <nlohmann/json.hpp>

#include "core/descriptors/Descriptor.hpp"
#include "core/kernels/Kernel.hpp"
#include "core/cutoff/CutoffFunction.hpp"
#include "data/kernels/TwoBodyKernelIndex.hpp"
#include "data/params/TwoBodyDescriptorParams.hpp"
#include "memory/MatrixBlock.hpp"

namespace jgap {

    class TwoBodyDescriptor : public Descriptor {
    public:
        explicit TwoBodyDescriptor(TwoBodyDescriptorParams params);

        ~TwoBodyDescriptor() override = default;

        double getCutoff() override {return _params.cutoff;};

        size_t nSparsePoints() override;

        nlohmann::json serialize() override;

        void setSparsePoints(const vector<AtomicStructure> &fromData) override;

        void setSparsePoints(map<SpeciesPair, vector<double>> sparsePointsPerSpeciesPair) {
            _sparsePointsPerSpeciesPair = std::move(sparsePointsPerSpeciesPair);
        }

        vector<Covariance> covariate(const AtomicStructure &atomicStructure) override;

        vector<pair<size_t, shared_ptr<MatrixBlock>>> selfCovariate() override;

    private:
        TwoBodyDescriptorParams _params;
        shared_ptr<Kernel<double, TwoBodyKernelIndex>> _kernel;
        shared_ptr<CutoffFunction> _cutoffFunction;

        map<SpeciesPair, vector<double>> _sparsePointsPerSpeciesPair;
        string _name;

        void sparsifyTrueUniform(const map<SpeciesPair, vector<double>> &all2b);
        void sparsifyQuipUniform(const map<SpeciesPair, vector<double>>& all2b);

        [[nodiscard]]
        map<SpeciesPair, TwoBodyKernelIndex> doIndex(const AtomicStructure &atomicStructure) const;
    };
}

#endif //TWOBODYDESCRIPTOR_HPP
