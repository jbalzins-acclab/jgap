#ifndef JGAP_ATOMSWITHNEIGHBOURLIST_HPP
#define JGAP_ATOMSWITHNEIGHBOURLIST_HPP
#include "NeighbourList.hpp"
#include "core/atomic/Atoms.hpp"

namespace jgap {
    struct AtomsWithNeighbourList {
        Atoms atoms;
        NeighbourList neighbour_list;
    };
}

#endif
