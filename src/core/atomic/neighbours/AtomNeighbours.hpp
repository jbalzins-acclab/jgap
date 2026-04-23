#ifndef JGAP_ATOMNEIGHBOURS_HPP
#define JGAP_ATOMNEIGHBOURS_HPP

#include "NeighbourData.hpp"
#include "core/atomic/geometry/SeparationTable.hpp"

namespace jgap {
    class AtomNeighboursData {


        std::map<Species, NeighbourData> neighbours;
        std::map<Species, SeparationTable> separations_between_neighbours;
        Real cutoff;
    };
}

#endif