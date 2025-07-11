#ifndef SAMPLESPACE3BSPARSIFIER_HPP
#define SAMPLESPACE3BSPARSIFIER_HPP

#include <nlohmann/json.hpp>

#include "core/descriptors/sparsification/Sparsifier.hpp"
#include "data/BasicDataTypes.hpp"

namespace jgap {

    using PointsPerSpeciesTriplet = map<SpeciesTriplet, vector<Vector3>>;

    class PerSpecies3bSparsifier : public Sparsifier<PointsPerSpeciesTriplet, vector<AtomicStructure>> {
    public:
        explicit PerSpecies3bSparsifier(nlohmann::json params);
        explicit PerSpecies3bSparsifier(const optional<vector<SpeciesTriplet>> &speciesTriplets)
            : _speciesTriplets(speciesTriplets) {}
        ~PerSpecies3bSparsifier() override = default;
        PointsPerSpeciesTriplet sparsifyFromData(const vector<AtomicStructure> &fromData) override;

    protected:
        optional<vector<SpeciesTriplet>> _speciesTriplets;

        virtual PointsPerSpeciesTriplet sparsifyFromTripletsInData(PointsPerSpeciesTriplet& tripletsInData) = 0;
    };
}

#endif //SAMPLESPACE3bSPARSIFIER_HPP
