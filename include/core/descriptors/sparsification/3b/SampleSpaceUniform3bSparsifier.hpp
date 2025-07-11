#ifndef SAMPLESPACEUNIFORM3BSPARSIFIER_HPP
#define SAMPLESPACEUNIFORM3BSPARSIFIER_HPP

#include "PerSpecies3bSparsifier.hpp"
#include "io/parse/ParserRegistry.hpp"

namespace jgap {
    class SampleSpaceUniform3bSparsifier : public PerSpecies3bSparsifier {
    public:
        explicit SampleSpaceUniform3bSparsifier(const nlohmann::json& params);

    protected:
        PointsPerSpeciesTriplet sparsifyFromTripletsInData(PointsPerSpeciesTriplet& all3b) override;

    private:
        size_t _nSparsePointsPerTriplet;
        optional<array<size_t, 3>> _nGridPoints;
        optional<array<array<double, 2/*from, to*/>, 3/*spec in each direction*/>> _sparseRanges;
    };

    REGISTER_PARSER("sample_space_uniform", PerSpecies3bSparsifier, SampleSpaceUniform3bSparsifier)
}

#endif //SAMPLESPACEUNIFORM3bSPARSIFIER_HPP
