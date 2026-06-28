#ifndef JGAP_TABULATIONDATA_HPP
#define JGAP_TABULATIONDATA_HPP

#include <map>

#include "ManyBodyGrids.hpp"
#include "TabulationParams.hpp"
#include "core/atomic/geometry/Cluster.hpp"
#include "core/atomic/species/SpeciesSet.hpp"
#include "core/splines/Grid.hpp"
#include "io/log/CurrentLogger.hpp"

namespace jgap {

    struct TabGapData {
        TabulationParams params;

        std::map<Species, Real> isolated_energies{};

        NBodyGrids<2, FullSymmetry> two_body_grids;
        NBodyGrids<3, NodeSymmetric> three_body_grids;

        std::vector<ManyBodyGrids<2>> eam_grids_vec{};

        TabGapData(const TabulationParams& params)
            : params(params),
              two_body_grids({0.0},
                             {params.max_cutoffs.forDim(2) / static_cast<Real>(params.n_grid_2b)},
                             {params.n_grid_2b}),
              three_body_grids({params.r_min_3b, params.r_min_3b, -1.0},
                               {
                                   params.max_cutoffs.forDim(3) / static_cast<Real>(params.n_grid_3b[0]),
                                    params.max_cutoffs.forDim(3) / static_cast<Real>(params.n_grid_3b[1]),
                                    2.0 / static_cast<Real>(params.n_grid_3b[2])
                               },
                               params.n_grid_3b) {
        }

        ManyBodyGrids<2>& newEamGrid(const Species& central_atom_species) {

            NBodyGrids<2, NodeSymmetric> aggregator_grids(
                {0.0},
                {params.max_cutoffs.per_cluster_size.at(2) / static_cast<Real>(params.n_grid_2b)},
                {params.n_grid_2b}
            );

            eam_grids_vec.emplace_back(
                central_atom_species,
                aggregator_grids,
                Grid<1>{
                    std::array{params.n_grid_2b},
                    std::array{params.max_eam_density / static_cast<Real>(params.n_grid_2b)},
                    std::array{0.0},
                }
            );
            return eam_grids_vec.back();
        }
    };

    struct TabulationData : TabGapData {

        template<size_t N>
        static Cluster<N> gridPosAsCluster(std::array<Real, Cluster<N>::NSeparations> grid_pos);

        explicit TabulationData(const TabulationParams& params) : TabGapData(params) {}
    };

    template<size_t N>
    Cluster<N> TabulationData::gridPosAsCluster(std::array<Real, Cluster<N>::NSeparations> grid_pos) {
        if constexpr (N == 2) {
            Cluster<2> res{};
            res.separationBetween(0, 1) = grid_pos[0];
            return res;
        }
        if constexpr (N == 3) {
            Cluster<3> res{};

            Real r01 = grid_pos[0];
            Real r02 = grid_pos[1];

            Real cos12 = grid_pos[2];
            Real r12 = sqrt(r01 * r01 + r02 * r02 - 2 * r01 * r02 * cos12);

            res.separationBetween(0, 1) = r01;
            res.separationBetween(0, 2) = r02;
            res.separationBetween(1, 2) = r12;
            return res;
        }

        // not static error to prevent in-progress N>3 descriptors from not compiling
        JGAP_LOG_AND_THROW("{}-Body tabulation not implemented");
    }
}

#endif