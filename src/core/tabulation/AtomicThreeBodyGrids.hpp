#ifndef JGAP_ATOMICTHREEBODYGRIDS_HPP
#define JGAP_ATOMICTHREEBODYGRIDS_HPP

#include <map>
#include "core/atomic/species/composition/Species3AtomicSorted.hpp"
#include "core/splines/Grid.hpp"

namespace jgap {
    template<size_t ValueDim = 1>
    struct AtomicThreeBodyGrids {
        static constexpr size_t DegreesOfFreedom = 3;
        static constexpr size_t GridDim = DegreesOfFreedom + (ValueDim > 1 ? 1 : 0);

        std::map<Species3AtomicSorted, Grid<GridDim>> value_grids;

        AtomicThreeBodyGrids(const std::array<Real, DegreesOfFreedom>& origin,
                             const std::array<Real, DegreesOfFreedom>& spacing,
                             const std::array<size_t, DegreesOfFreedom>& dims)
            : origin(origin), spacing(spacing), dims(dims) {}

        Grid<GridDim>& getValueGrid(const Species3AtomicSorted& species_set) {
            if (isPlaceholder()) {
                JGAP_LOG_AND_THROW("Spacing is zero in one of the dimensions,"
                                   "check that the max cutoff is non-zero");
            }

            auto it = value_grids.find(species_set);
            if (it == value_grids.end()) {
                std::array<size_t, GridDim> grid_dims;
                std::array<Real, GridDim> grid_spacing;
                std::array<Real, GridDim> grid_origin;
                for (size_t i = 0; i < DegreesOfFreedom; ++i) {
                    grid_dims[i] = dims[i];
                    grid_spacing[i] = spacing[i];
                    grid_origin[i] = origin[i];
                }
                if constexpr (ValueDim > 1) {
                    grid_dims[DegreesOfFreedom] = ValueDim;
                    grid_spacing[DegreesOfFreedom] = 1.0;
                    grid_origin[DegreesOfFreedom] = 0.0;
                }

                auto [new_it, inserted] = value_grids.emplace(
                    species_set,
                    Grid<GridDim>(grid_dims, grid_spacing, grid_origin)
                );
                return new_it->second;
            }
            return it->second;
        }

        bool hasGridFor(const Species3AtomicSorted& species_set) const {
            return value_grids.contains(species_set);
        }

        bool isPlaceholder() const {
            bool res = false;
            for (size_t dim{}; dim < DegreesOfFreedom; dim++) {
                res |= (dims[dim] == 0);
            }
            return res;
        }

    private:
        std::array<Real, DegreesOfFreedom> origin;
        std::array<Real, DegreesOfFreedom> spacing;
        std::array<size_t, DegreesOfFreedom> dims;
    };
}

#endif
