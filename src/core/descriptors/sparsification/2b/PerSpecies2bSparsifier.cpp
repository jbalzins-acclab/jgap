#include "core/descriptors/sparsification/2b/PerSpecies2bSparsifier.hpp"

#include "utils/Utils.hpp"

namespace jgap {
    PerSpecies2bSparsifier::PerSpecies2bSparsifier(nlohmann::json params) {
        if (params.contains("species")) {
            vector<SpeciesPair> sps;
            for (const auto &speciesPair : params["species"]) {
                sps.push_back(SpeciesPair{split(speciesPair, ',')[0], split(speciesPair, ',')[1]});
            }
            _speciesPairs = sps;
        }
    }

    PointsPerSpeciesPair PerSpecies2bSparsifier::sparsifyFromData(const vector<AtomicStructure> &fromData) {

        map<SpeciesPair, vector<double>> all2b;
        // Explicitly defined pairs
        if (_speciesPairs.has_value()) {
            for (auto &speciesPair : _speciesPairs.value()) {
                all2b[speciesPair] = {};
            }
        }

        for (const auto& structure: fromData) {
            for (const auto& atomData : structure.atoms) {
                if (!atomData.neighbours.has_value()) {
                    CurrentLogger::get() -> error("Neighbour list missing | 2b", true);
                }

                for (const NeighbourData& neighbour : atomData.neighbours.value()) {

                    auto species = SpeciesPair(atomData.species, structure.atoms[neighbour.index].species);
                    if (!all2b.contains(species)) {
                        if (_speciesPairs.has_value()) continue;
                        // pairs-not explicitly defined => from data
                        all2b[species] = {};
                    }
                    if (all2b.contains(species)) {
                        all2b[species].push_back(neighbour.distance);
                    }
                }
            }
        }

        return sparsifyFromDistancesInData(all2b);
    }
}
