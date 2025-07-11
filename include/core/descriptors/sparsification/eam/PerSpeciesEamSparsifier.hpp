#ifndef PERSPECIESEAMBSPARSIFIER_HPP
#define PERSPECIESEAMBSPARSIFIER_HPP

#include <nlohmann/json.hpp>

#include "core/descriptors/sparsification/Sparsifier.hpp"
#include "data/BasicDataTypes.hpp"
#include "data/kernels/EamKernelIndex.hpp"

namespace jgap {

    using DensitiesPerSpecies = map<Species, vector<double>>;

    class PerSpeciesEamSparsifier : public Sparsifier<DensitiesPerSpecies, vector<EamKernelIndex>> {
    public:
        explicit PerSpeciesEamSparsifier(nlohmann::json params);
        explicit PerSpeciesEamSparsifier(const optional<vector<Species>> &species) : _species(species) {}

        DensitiesPerSpecies sparsifyFromData(const vector<EamKernelIndex> &fromData) override;

    protected:
        optional<vector<Species>> _species;

        virtual DensitiesPerSpecies sparsifyFromDensitiesInData(DensitiesPerSpecies& densitiesInData) = 0;
    };
}


#endif