#ifndef JGAP_MANYBODYGRIDS2_HPP
#define JGAP_MANYBODYGRIDS2_HPP
#include "AtomicTwoBodyGrids.hpp"
#include "core/splines/Grid.hpp"

namespace jgap {
    template<size_t AggregatorDim = 1, size_t ValueDim = 1>
    struct ManyBodyGrids2 {

        Species central_atom_species;

        AtomicTwoBodyGrids<AggregatorDim> aggregator_grids;
        Grid<AggregatorDim + (ValueDim > 1 ? 1 : 0)> value_grid;
    };
}

#endif
