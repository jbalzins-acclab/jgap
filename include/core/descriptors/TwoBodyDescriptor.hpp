#ifndef TWOBODYDESCRIPTOR_HPP
#define TWOBODYDESCRIPTOR_HPP

#include <ranges>
#include <utility>
#include <nlohmann/json.hpp>

#include "core/descriptors/Descriptor.hpp"
#include "core/descriptors/kernels/Kernel.hpp"
#include "core/cutoff/CutoffFunction.hpp"
#include "data/kernels/TwoBodyKernelIndex.hpp"
#include "io/parse/ParserRegistry.hpp"
#include "core/descriptors/kernels/TwoBodySE.hpp"
#include "memory/MatrixBlock.hpp"
#include "core/descriptors/sparsification/2b/PerSpecies2bSparsifier.hpp"

namespace jgap {

    class TwoBodyDescriptor : public Descriptor {
    public:

        ~TwoBodyDescriptor() override = default;

        explicit TwoBodyDescriptor(const nlohmann::json& params);
        nlohmann::json serialize() override;
        string getType() override { return "2b"; };

        double getCutoff() override { return _cutoff; };

        size_t nSparsePoints() override;
        void setSparsePoints(const vector<AtomicStructure> &fromData) override;
        void setSparsePoints(map<SpeciesPair, vector<double>> sparsePointsPerSpeciesPair) {
            _sparsePointsPerSpeciesPair = std::move(sparsePointsPerSpeciesPair);
        }

        vector<Covariance> covariate(const AtomicStructure &atomicStructure) override;
        vector<pair<size_t, shared_ptr<MatrixBlock>>> selfCovariate() override;

        TabulationData tabulate(const TabulationParams &params) override;

    private:
        double _cutoff;
        shared_ptr<TwoBodyKernel> _kernel;
        shared_ptr<PerSpecies2bSparsifier> _sparsifier;

        map<SpeciesPair, vector<double>> _sparsePointsPerSpeciesPair;

        [[nodiscard]]
        map<SpeciesPair, TwoBodyKernelIndex> doIndex(const AtomicStructure &atomicStructure) const;
    };

    REGISTER_PARSER("2b", Descriptor, TwoBodyDescriptor)
}

#endif //TWOBODYDESCRIPTOR_HPP
