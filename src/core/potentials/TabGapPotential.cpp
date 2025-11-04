#include "core/potentials/TabGapPotential.hpp"

#include "data/AtomicStructure.hpp"
#include "data/PredictionData.hpp"

namespace jgap {
    TabGapPotential::TabGapPotential(const nlohmann::json &params) {
    }

    Predictions TabGapPotential::predict(const AtomicStructure &structure) {

        double energy = 0.0;
        vector forces(structure.size(), Vector3(0.0, 0.0, 0.0));
        array<Vector3, 3> virials{};

        for (const auto& atom: structure) {
            energy += _isolatedEnergies[atom.species()];

            vector<size_t> ;
            for (const NeighbourData& neighbour: atom.neighbours().size()) {

            }
        }

        return {energy, forces, virials};
    }

    void TabGapPotential::tabulate(TabulationData &table) {
    }

    void TabGapPotential::saveFiles(string) {
    }
}
