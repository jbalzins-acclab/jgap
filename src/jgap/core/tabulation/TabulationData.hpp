#ifndef JGAP_TABULATIONDATA_HPP
#define JGAP_TABULATIONDATA_HPP

#include <map>
#include <optional>
#include <utility>

#include "AtomicThreeBodyGrids.hpp"
#include "AtomicTwoBodyGrids.hpp"
#include "ManyBodyGrids2.hpp"
#include "TabulationParams.hpp"
#include "TwoBodyGrids.hpp"
#include "jgap/core/atomic/geometry/Cluster2.hpp"
#include "jgap/core/atomic/geometry/Cluster3.hpp"
#include "jgap/core/atomic/species/composition/Species3AtomicSorted.hpp"
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
                {0.0_r}, {params.max_cutoffs.forDim(2) / static_cast<Real>(params.n_grid_2b - 1)}, {params.n_grid_2b}
            ),
            three_body_grids(
                {params.r_min_3b, params.r_min_3b, -1.0_r},
                {(params.max_cutoffs.forDim(3) - params.r_min_3b) / static_cast<Real>(params.n_grid_3b[0] - 1),
                 (params.max_cutoffs.forDim(3) - params.r_min_3b) / static_cast<Real>(params.n_grid_3b[1] - 1),
                 2.0_r / static_cast<Real>(params.n_grid_3b[2] - 1)},
                params.n_grid_3b
            ) {}

        ManyBodyGrids2<1, 1>& newEamGrid(const Species& central_atom_species) {
            AtomicTwoBodyGrids<1> aggregator_grids(
                {0.0_r},
                {params.max_cutoffs.per_cluster_size.at(2) / static_cast<Real>(params.n_grid_2b - 1)},
                {params.n_grid_2b}
            );

            eam_grids_vec.emplace_back(
                central_atom_species,
                aggregator_grids,
                Grid<1>{
                    std::array{params.n_grid_2b},
                    std::array{params.max_eam_density / static_cast<Real>(params.n_grid_2b - 1)},
                    std::array{0.0_r},
                }
            );
            return eam_grids_vec.back();
        }
    };

    struct TabulationData : TabGapData {
        static Cluster2 gridPosAsCluster2(std::array<Real, 1> grid_pos);
        static std::pair<Cluster3, std::optional<Cluster3>> gridPosAsCluster3(
            std::array<Real, 3> grid_pos, const Species3AtomicSorted& species
        );

        explicit TabulationData(const TabulationParams& params) : TabGapData(params) {}

    private:
        static Cluster3 gridPosAsCluster3(std::array<Real, 3> grid_pos);
    };

    inline Cluster2 TabulationData::gridPosAsCluster2(std::array<Real, 1> grid_pos) {
        Cluster2 res{};
        res.separation01.magnitude = grid_pos[0];
        return res;
    }

    inline Cluster3 TabulationData::gridPosAsCluster3(std::array<Real, 3> grid_pos) {
        Cluster3 res{};

        Real r01 = grid_pos[0];
        Real r02 = grid_pos[1];

        Real cos12 = std::clamp(grid_pos[2], -1.0_r, 1.0_r);
        Real term = r01 * r01 + r02 * r02 - 2 * r01 * r02 * cos12;
        Real r12 = sqrt(std::max(0.0_r, term));

        res.separation01.magnitude = r01;
        res.separation02.magnitude = r02;
        res.separation12.magnitude = r12;
        return res;
    }

    inline std::pair<Cluster3, std::optional<Cluster3>> TabulationData::gridPosAsCluster3(
        std::array<Real, 3> grid_pos, const Species3AtomicSorted& species
    ) {
        Cluster3 cluster1 = gridPosAsCluster3(grid_pos);
        if (species.nodes[0] == species.nodes[1]) {
            Cluster3 cluster2 = cluster1;
            cluster2.separation01.magnitude = cluster1.separation02.magnitude;
            cluster2.separation02.magnitude = cluster1.separation01.magnitude;
            return {cluster1, cluster2};
        }
        return {cluster1, std::nullopt};
    }
}

#endif
