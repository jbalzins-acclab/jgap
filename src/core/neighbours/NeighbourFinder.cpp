#include "core/neighbours/NeighbourFinder.hpp"

#include <tbb/parallel_for_each.h>
#include <iostream>
using namespace std;

namespace jgap {
    void NeighbourFinder::findNeighbours(vector<AtomicStructure> &structures, double cutoff) {
        tbb::parallel_for_each(
            structures.begin(), structures.end(),
            [cutoff](AtomicStructure &structure) { findNeighbours(structure, cutoff); }
        );
    }

    tuple<int, int, int> findMaxRep(const AtomicStructure& structure, const double cutoff) {
        const Vector3 side1 = structure.lattice[0],
                      side2 = structure.lattice[1],
                      side3 = structure.lattice[2];

        tuple<int, int, int> maxRep = {
            cutoff / side1.norm() + 2,
            cutoff / side2.norm() + 2,
            cutoff / side3.norm() + 2
        };

        // triclinic
        if (side1.dot(side2) > 1e-6 || side1.dot(side3) > 1e-6 || side2.dot(side3) > 1e-6) {
            get<0>(maxRep) = static_cast<int>(cutoff / side1.aproject(side2, side3)) + 2;
            get<1>(maxRep) = static_cast<int>(cutoff / side2.aproject(side1, side3)) + 2;
            get<2>(maxRep) = static_cast<int>(cutoff / side3.aproject(side2, side1)) + 2;
        }

        return maxRep;
    }

    void NeighbourFinder::findNeighbours(AtomicStructure& structure, const double cutoff) {

        /*
        vector<Vector3> corners;
        for (int mask = 0; mask < 8; mask++) {
            auto corner =
                structure->latticeVectors[0] * ((mask & 1) != 0)
                + structure->latticeVectors[1] * ((mask & 2) != 0)
                + structure->latticeVectors[2] * ((mask & 4) != 0);
            corners.push_back(corner);
        }
        */

        const auto maxRep = findMaxRep(structure, cutoff);

        vector<Vector3> possibleOffsets;
        const auto zeroVec = Vector3(0, 0, 0);

        for (int rep0 = -get<0>(maxRep); rep0 <= get<0>(maxRep); rep0++) {
            for (int rep1 = -get<1>(maxRep); rep1 <= get<1>(maxRep); rep1++) {
                for (int rep2 = -get<2>(maxRep); rep2 <= get<2>(maxRep); rep2++) {

                    auto offset = zeroVec
                        + structure.lattice[0] * rep0
                        + structure.lattice[1] * rep1
                        + structure.lattice[2] * rep2;

                    possibleOffsets.push_back(offset);
                }
            }
        }

        // avoid heap corruption
        if (!structure.neighbours.has_value()) {
            structure.neighbours = vector(structure.size(), NeighboursData());
        } else {
            for (size_t i = 0; i < structure.size(); i++) {
                (*structure.neighbours)[i].clear();
            }
        }

        for (size_t i = 0; i < structure.size(); i++) {
            for (size_t j = i; j < structure.size(); j++) {
                for (const auto &offset : possibleOffsets) {
                    auto dist = (structure[i].position() - (structure[j].position() + offset)).norm();
                    if (0 < dist && dist <= cutoff) {
                        (*structure.neighbours)[i].push_back({.index=j, .offset=offset, .distance = dist});
                        if (i != j) {
                            (*structure.neighbours)[j].push_back({.index=i, .offset=offset * -1, .distance = dist});
                        }
                    }
                }
            }
        }
    }
}
