#ifndef PERSPECIES2BSPARSIFIER_HPP
#define PERSPECIES2BSPARSIFIER_HPP

#include <nlohmann/json.hpp>

#include "core/descriptors/sparsification/Sparsifier.hpp"
#include "data/BasicDataTypes.hpp"

namespace jgap {

    using PointsPerSpeciesPair = map<SpeciesPair, vector<double>>;

    class PerSpecies2bSparsifier : public Sparsifier<PointsPerSpeciesPair, vector<AtomicStructure>> {
    public:
        explicit PerSpecies2bSparsifier(nlohmann::json params);
        explicit PerSpecies2bSparsifier(const optional<vector<SpeciesPair>> &speciesPairs)
            : _speciesPairs(speciesPairs) {}
        ~PerSpecies2bSparsifier() override = default;
        PointsPerSpeciesPair sparsifyFromData(const vector<AtomicStructure> &fromData) override;

    protected:
        optional<vector<SpeciesPair>> _speciesPairs;

        virtual PointsPerSpeciesPair sparsifyFromDistancesInData(PointsPerSpeciesPair& distancesInData) = 0;
    };
}

#endif
