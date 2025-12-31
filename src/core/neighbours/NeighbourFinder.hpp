#ifndef JGAP_NEIGHBOURFINDER_HPP
#define JGAP_NEIGHBOURFINDER_HPP

#include "data/Vector3.hpp"

#include <vector>

#include "../../data/atomic/AtomicStructure.hpp"

namespace jgap {
    class NeighbourFinder {
    public:
        static void findNeighbours(std::vector<AtomicStructure> &structures, double cutoff);
        static void findNeighbours(AtomicStructure& structure, double cutoff);
    };
}

#endif
