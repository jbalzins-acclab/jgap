#ifndef JGAP_MANYBODYGRIDS_HPP
#define JGAP_MANYBODYGRIDS_HPP
#include "NBodyGrids.hpp"
#include "core/splines/Grid.hpp"

namespace jgap {
    template<size_t NBodies, size_t AggregatorDim = 1, size_t ValueDim = 1>
    struct ManyBodyGrids {

        Species central_atom_species;

        NBodyGrids<NBodies, NodeSymmetric, AggregatorDim> aggregator_grids;
        Grid<AggregatorDim> value_grid;
    };
}

#endif
