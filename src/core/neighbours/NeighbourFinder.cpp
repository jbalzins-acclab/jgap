#include "core/neighbours/NeighbourFinder.hpp"

#include <tbb/parallel_for_each.h>
#include <iostream>
#include <algorithm>

namespace jgap {
    void NeighbourFinder::findNeighbours(std::vector<AtomicStructure> &structures, double cutoff) {
        tbb::parallel_for_each(
            structures.begin(), structures.end(),
            [cutoff](AtomicStructure &structure) { findNeighbours(structure, cutoff); }
        );
    }

    std::array<int, 3> findMaxRep(const AtomicStructure& structure, const double cutoff) {
        const Vector3 side1 = structure.lattice_vectors[0],
                      side2 = structure.lattice_vectors[1],
                      side3 = structure.lattice_vectors[2];

        std::array maxRep = {
            static_cast<int>(cutoff / side1.len() + 2),
            static_cast<int>(cutoff / side2.len() + 2),
            static_cast<int>(cutoff / side3.len() + 2)
        };

        // triclinic
        if (side1.dot(side2) > 1e-6 || side1.dot(side3) > 1e-6 || side2.dot(side3) > 1e-6) {
            maxRep[0] = static_cast<int>(cutoff / side1.aproject(side2, side3)) + 2;
            maxRep[1] = static_cast<int>(cutoff / side2.aproject(side1, side3)) + 2;
            maxRep[2] = static_cast<int>(cutoff / side3.aproject(side2, side1)) + 2;
        }

        return maxRep;
    }

    void NeighbourFinder::findNeighbours(AtomicStructure& structure, const double cutoff) {

        const auto maxRep = findMaxRep(structure, cutoff);

        std::vector<Vector3> possibleOffsets;
        const auto zeroVec = Vector3(0, 0, 0);

        for (int rep0 = -get<0>(maxRep); rep0 <= get<0>(maxRep); rep0++) {
            for (int rep1 = -get<1>(maxRep); rep1 <= get<1>(maxRep); rep1++) {
                for (int rep2 = -get<2>(maxRep); rep2 <= get<2>(maxRep); rep2++) {

                    auto offset = zeroVec
                        + structure.lattice_vectors[0] * rep0
                        + structure.lattice_vectors[1] * rep1
                        + structure.lattice_vectors[2] * rep2;

                    possibleOffsets.push_back(offset);
                }
            }
        }

        // avoid heap corruption
        if (!structure.neighbours_ascending_separation.has_value()) {
            structure.neighbours_ascending_separation = std::vector(structure.size(), NeighboursData());
        } else {
            for (size_t i = 0; i < structure.size(); i++) {
                (*structure.neighbours_ascending_separation)[i].clear();
            }
        }

        for (size_t i = 0; i < structure.size(); i++) {
            for (size_t j = i; j < structure.size(); j++) {
                for (const auto &offset : possibleOffsets) {
                    auto dist = (structure[i].position() - (structure[j].position() + offset)).len();
                    if (0 < dist && dist <= cutoff) {
                        (*structure.neighbours_ascending_separation)[i].push_back({
                            .index=j, .offset=offset, .distance = dist
                        });
                        if (i != j) {
                            (*structure.neighbours_ascending_separation)[j].push_back({
                                .index = i, .offset = offset * -1, .distance = dist
                            });
                        }
                    }
                }
            }
        }
        
        for (auto& neighbours_data : *structure.neighbours_ascending_separation) {
            std::ranges::sort(neighbours_data, [](auto& a, auto& b) { return a.distance < b.distance; });
        }
    }
}
