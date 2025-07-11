#include "core/descriptors/sparsification/3b/PerSpecies3bSparsifier.hpp"

namespace jgap {
    PerSpecies3bSparsifier::PerSpecies3bSparsifier(nlohmann::json params) {
        if (params.contains("species")) {
            _speciesTriplets = vector<SpeciesTriplet>();
            for (const auto &triplet: params["species"]) {
                _speciesTriplets->push_back(
                    SpeciesTriplet{
                        triplet["root"].get<string>(),
                        {triplet["node1"].get<string>(), triplet["node2"].get<string>()}
                    }
                );
            }
        }
    }

    PointsPerSpeciesTriplet PerSpecies3bSparsifier::sparsifyFromData(const vector<AtomicStructure> &fromData) {

        map<SpeciesTriplet, vector<Vector3>> all3b;

        if (_speciesTriplets.has_value()) {
            for (auto &speciesTriplet : _speciesTriplets.value()) {
                all3b[speciesTriplet] = {};
            }
        }

        for (const auto& structure: fromData) {
            for (const auto& atom0 : structure.atoms) {
                if (!atom0.neighbours.has_value()) {
                    CurrentLogger::get() -> error("Neighbour list missing | 3b");
                    throw runtime_error("Neighbour list missing | 3b");
                }

                for (size_t idx1 = 0; idx1 < atom0.neighbours->size(); idx1++) {

                    auto neighbour1 = atom0.neighbours->at(idx1);
                    // if (neighbour1.distance > _cutoff) continue;

                    for (size_t idx2 = idx1+1; idx2 < atom0.neighbours->size(); idx2++) {
                        auto neighbour2 = atom0.neighbours->at(idx2);

                        // if (neighbour2.distance > _cutoff) continue;
                        if (neighbour1.index == neighbour2.index && neighbour1.offset == neighbour2.offset) continue; //??

                        auto atom1 = structure.atoms[neighbour1.index], atom2 = structure.atoms[neighbour2.index];

                        auto species = SpeciesTriplet{
                            atom0.species,
                            {atom1.species, atom2.species}
                        };

                        if (!all3b.contains(species)) {
                            if (_speciesTriplets.has_value()) continue;
                            // pairs-not explicitly defined => from data
                            all3b[species] = {};
                        }
                        if (all3b.contains(species)) {
                            all3b[species].emplace_back(
                                neighbour1.distance,
                                neighbour2.distance,
                                (atom1.position + neighbour1.offset - (atom2.position + neighbour2.offset)).norm()
                            );
                        }
                    }
                }
            }
        }

        return sparsifyFromTripletsInData(all3b);
    }
}