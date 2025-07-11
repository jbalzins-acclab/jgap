#ifndef SAMPLESPACEUNIFORMEAMSPARSIFIER_HPP
#define SAMPLESPACEUNIFORMEAMSPARSIFIER_HPP

#include "core/descriptors/sparsification/eam/PerSpeciesEamSparsifier.hpp"
#include "io/parse/ParserRegistry.hpp"

namespace jgap {
    class SampleSpaceUniformEamSparsifier : public PerSpeciesEamSparsifier {
    public:
        SampleSpaceUniformEamSparsifier(const nlohmann::json& params);

    protected:
        size_t _nSparsePoints;
        map<Species, array<double, 2>> _sparseRanges;

        DensitiesPerSpecies sparsifyFromDensitiesInData(DensitiesPerSpecies& densitiesInData) override;
    };

    REGISTER_PARSER("sample_space_uniform", PerSpeciesEamSparsifier, SampleSpaceUniformEamSparsifier)
}

#endif //SAMPLESPACEUNIFORMEAMSPARSIFIER_HPP
