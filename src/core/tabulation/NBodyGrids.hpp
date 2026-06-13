#ifndef JGAP_NBODYGRIDS_HPP
#define JGAP_NBODYGRIDS_HPP

#include "core/atomic/species/SpeciesSet.hpp"
#include "core/splines/Grid.hpp"
#include "core/atomic/geometry/Cluster.hpp"

namespace jgap {
    template<size_t N, ClusterSymmetry ClusterSym, size_t ValueDim = 1>
    requires(N > 1 && ValueDim > 0)
    struct NBodyGrids;

    template<size_t N, ClusterSymmetry ClusterSym>
    struct NBodyGrids<N, ClusterSym> {
        static constexpr size_t DegreesOfFreedom{Cluster<N>::NSeparations};

        std::map<SpeciesSet<N, ClusterSym>, Grid<DegreesOfFreedom>> value_grids;

        NBodyGrids(const std::array<Real, DegreesOfFreedom>& origin,
                   const std::array<Real, DegreesOfFreedom>& spacing,
                   const std::array<size_t, DegreesOfFreedom>& dims)
            : origin(origin), spacing(spacing), dims(dims) {}

        Grid<DegreesOfFreedom>& getValueGrid(const SpeciesSet<N, ClusterSym>& species_set) {

            if (isPlaceholder()) {
                JGAP_LOG_AND_THROW("Spacing is zero in one of the {} dimensions,"
                                   "check that the max cutoff is non-zero",
                                   DegreesOfFreedom);
            }

            auto it = value_grids.find(species_set);
            if (it == value_grids.end()) {
                auto [new_it, inserted] = value_grids.emplace(
                    species_set,
                    Grid<DegreesOfFreedom>(dims, spacing, origin)
                );
                return new_it->second;
            }
            return it->second;
        }

        bool hasGridFor(const SpeciesSet<N, ClusterSym>& species_set) const {
            return value_grids.contains(species_set);
        }

        bool isPlaceholder() const {
            bool res = false;
            for (size_t dim{}; dim < DegreesOfFreedom; dim++) {
                res |= (dims[dim] == 0.0);
            }
            return res;
        }

    private:
        std::array<Real, DegreesOfFreedom> origin;
        const std::array<Real, DegreesOfFreedom> spacing;
        const std::array<size_t, DegreesOfFreedom> dims;
    };
}

#endif
