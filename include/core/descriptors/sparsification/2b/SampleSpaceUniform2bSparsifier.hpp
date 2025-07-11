#ifndef SAMPLESPACEUNIFORM2BSPARSIFIER_HPP
#define SAMPLESPACEUNIFORM2BSPARSIFIER_HPP

#include "PerSpecies2bSparsifier.hpp"
#include "io/parse/ParserRegistry.hpp"

namespace jgap {
    class SampleSpaceUniform2bSparsifier : public PerSpecies2bSparsifier {
    public:
        explicit SampleSpaceUniform2bSparsifier(const nlohmann::json &params);

    protected:
        PointsPerSpeciesPair sparsifyFromDistancesInData(PointsPerSpeciesPair& distancesInData) override;

    private:
        size_t _nSparsePointsPerSpeciesPair;
        optional<array<double, 2>> _sparseRange;
    };

    REGISTER_PARSER("sample_space_uniform", PerSpecies2bSparsifier, SampleSpaceUniform2bSparsifier)
}

#endif //SAMPLESPACEUNIFORM2BSPARSIFIER_HPP
