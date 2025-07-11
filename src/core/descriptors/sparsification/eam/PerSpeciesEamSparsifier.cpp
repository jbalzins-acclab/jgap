#include "core/descriptors/sparsification/eam/PerSpeciesEamSparsifier.hpp"

namespace jgap {
    PerSpeciesEamSparsifier::PerSpeciesEamSparsifier(nlohmann::json params) {
        _species = {};
        if (params.contains("species")) {
            vector<Species> speciesArr;
            for (const auto& species : params["species"]) {
                speciesArr.push_back(species.get<string>());
            }
            _species = speciesArr;
        }
    }

    DensitiesPerSpecies PerSpeciesEamSparsifier::sparsifyFromData(const vector<EamKernelIndex> &fromData) {

        DensitiesPerSpecies densities{};

        for (const auto &index : fromData) {
            for (const auto &[species, densityDatas]: index) {
                for (const auto &densityData : densityDatas) {
                    densities[species].push_back(densityData.density);
                }
            }
        }

        return sparsifyFromDensitiesInData(densities);
    }
}
