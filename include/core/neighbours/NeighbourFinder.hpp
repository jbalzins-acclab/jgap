#ifndef NEIGHBOURFINDER_HPP
#define NEIGHBOURFINDER_HPP

#include "data/Vector3.hpp"

#include <vector>
#include <memory>


namespace jgap {
    class NeighbourFinder {
    public:
        static void findNeighbours(vector<AtomicStructure> &structures, double cutoff);
        static void findNeighbours(AtomicStructure& structure, double cutoff);
    };
}

#endif //NEIGHBOURFINDER_HPP
