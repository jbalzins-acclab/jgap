#ifndef EAMDESCRIPTOR_HPP
#define EAMDESCRIPTOR_HPP

#include <nlohmann/json.hpp>

#include "core/descriptors/Descriptor.hpp"
#include "core/descriptors/kernels/Kernel.hpp"
#include "data/kernels/EamKernelIndex.hpp"
#include "eam/pair_functions/EamPairFunction.hpp"
#include "io/parse/ParserRegistry.hpp"
#include "kernels/EamSE.hpp"
#include "memory/MatrixBlock.hpp"
#include "core/descriptors/sparsification/eam/PerSpeciesEamSparsifier.hpp"

namespace jgap {

    class EamDescriptor : public Descriptor {
    public:
        explicit EamDescriptor(const nlohmann::json &params);
        nlohmann::json serialize() override;
        string getType() override { return "eam"; };

        ~EamDescriptor() override = default;

        void setSparsePoints(const vector<AtomicStructure> &fromData) override;

        void setSparsePoints(map<Species, vector<double>> sparsePoints) {
            _sparsePointsPerSpecies = std::move(sparsePoints);
        }

        size_t nSparsePoints() override;
        map<Species, vector<double>> getSparsePoints() { return _sparsePointsPerSpecies; }

        double getCutoff() override { return _cutoff; }

        vector<Covariance> covariate(const AtomicStructure &atomicStructure) override;
        vector<pair<size_t, shared_ptr<MatrixBlock>>> selfCovariate() override;

        TabulationData tabulate(const TabulationParams &params) override;

    private:
        double _cutoff;
        shared_ptr<EamKernel> _kernel;
        shared_ptr<PerSpeciesEamSparsifier> _sparsifier;

        shared_ptr<EamPairFunction> _defaultPairFunction;
        map<OrderedSpeciesPair, shared_ptr<EamPairFunction>> _pairFunctions;

        map<Species, vector<double>> _sparsePointsPerSpecies;

        [[nodiscard]] EamKernelIndex doIndex(const AtomicStructure &structure) const;
    };

    REGISTER_PARSER("eam", Descriptor, EamDescriptor)
}



#endif //EAMDESCRIPTOR_HPP
