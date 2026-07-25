#ifndef JGAP_TABULATIONDATA_HPP
#define JGAP_TABULATIONDATA_HPP

#include <map>

#include "AtomicThreeBodyGrids.hpp"
#include "AtomicTwoBodyGrids.hpp"
#include "ManyBodyGrids2.hpp"
#include "TabulationParams.hpp"
#include "TwoBodyGrids.hpp"
#include "jgap/core/atomic/geometry/Cluster2.hpp"
#include "jgap/core/atomic/geometry/Cluster3.hpp"
#include "jgap/core/splines/Grid.hpp"

namespace jgap {

    struct TabGapData {
        TabulationParams params;

        std::map<Species, Real> isolated_energies{};

        TwoBodyGrids<1> two_body_grids;
        AtomicThreeBodyGrids<1> three_body_grids;

        std::vector<ManyBodyGrids2<1, 1>> eam_grids_vec{};

        TabGapData(const TabulationParams& params) :
            params(params),
            two_body_grids(
                {0.0_r}, {params.max_cutoffs.forDim(2) / static_cast<Real>(params.n_grid_2b)}, {params.n_grid_2b}
            ),
            three_body_grids(
                {params.r_min_3b, params.r_min_3b, -1.0_r},
                {params.max_cutoffs.forDim(3) / static_cast<Real>(params.n_grid_3b[0]),
                 params.max_cutoffs.forDim(3) / static_cast<Real>(params.n_grid_3b[1]),
                 2.0_r / static_cast<Real>(params.n_grid_3b[2])},
                params.n_grid_3b
            ) {}

        ManyBodyGrids2<1, 1>& newEamGrid(const Species& central_atom_species) {
            AtomicTwoBodyGrids<1> aggregator_grids(
                {0.0_r}, {params.max_cutoffs.per_cluster_size.at(2) / static_cast<Real>(params.n_grid_2b)},
                {params.n_grid_2b}
            );

            eam_grids_vec.emplace_back(
                central_atom_species, aggregator_grids,
                Grid<1>{
                    std::array{params.n_grid_2b},
                    std::array{params.max_eam_density / static_cast<Real>(params.n_grid_2b)},
                    std::array{0.0_r},
                }
            );
            return eam_grids_vec.back();
        }
    };

    struct TabulationData : TabGapData {
        static Cluster2 gridPosAsCluster2(std::array<Real, 1> grid_pos);
        static Cluster3 gridPosAsCluster3(std::array<Real, 3> grid_pos);

        explicit TabulationData(const TabulationParams& params) : TabGapData(params) {}
    };

    inline Cluster2 TabulationData::gridPosAsCluster2(std::array<Real, 1> grid_pos) {
        Cluster2 res{};
        res.r01 = grid_pos[0];
        return res;
    }

    inline Cluster3 TabulationData::gridPosAsCluster3(std::array<Real, 3> grid_pos) {
        Cluster3 res{};

        Real r01 = grid_pos[0];
        Real r02 = grid_pos[1];

        Real cos12 = grid_pos[2];
        Real r12 = sqrt(r01 * r01 + r02 * r02 - 2 * r01 * r02 * cos12);

        res.separationBetween(0, 1) = r01;
        res.separationBetween(0, 2) = r02;
        res.separationBetween(1, 2) = r12;
        return res;
    }
}

#endif
