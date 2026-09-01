#include <cmath>
#include <gtest/gtest.h>

#include "jgap/core/atomic/Atoms.hpp"
#include "jgap/core/atomic/species/composition/Species2Sorted.hpp"
#include "jgap/core/atomic/species/composition/Species3AtomicSorted.hpp"
#include "jgap/core/cutoff/CosCutoff.hpp"
#include "jgap/core/kernels/SquaredExpKernel.hpp"
#include "jgap/core/potentials/gap/GapPotential.hpp"
#include "jgap/core/potentials/gap/component/ThreeBodyGapComponent.hpp"
#include "jgap/core/potentials/gap/component/ManyBodyGapComponent.hpp"
#include "jgap/core/potentials/gap/component/TwoBodyGapComponent.hpp"
#include "jgap/core/potentials/isolated/IsolatedAtomPotential.hpp"
#include "jgap/core/potentials/tabgap/TabGapPotential.hpp"
#include "jgap/core/transform/manybody/TwoBodySum.hpp"
#include "jgap/core/transform/nbody/2b/PairDistanceTransformation.hpp"
#include "jgap/core/transform/nbody/2b/eam/PolycutoffPairFunction.hpp"
#include "jgap/core/transform/nbody/3b/Angle3bTransformation.hpp"

using namespace jgap;

TEST(TestTabGapPotential, IsolatedAtomEnergy) {
    std::map<Species, Real> expected_iso{{"Fe", -4.2}, {"Ni", -3.1}};
    auto iso_pot = IsolatedAtomPotential(expected_iso);

    GapPotential gap;
    gap.optional_external_potential = ValuePtr<Potential>(iso_pot);

    TabulationData tables = gap.tabulate({.n_grid_2b = 100});
    TabGapPotential tabgap(tables);

    const auto& iso_energies = tabgap.getIsolatedAtomEnergies();
    EXPECT_EQ(iso_energies.size(), 2);
    EXPECT_DOUBLE_EQ(iso_energies.at(Species("Fe")), -4.2);
    EXPECT_DOUBLE_EQ(iso_energies.at(Species("Ni")), -3.1);

    Atoms atoms({{0, 0, 0}, {2.5, 0, 0}}, {Species("Fe"), Species("Ni")});
    auto energy_result = tabgap.calculateEnergy(atoms);
    EXPECT_NEAR(energy_result.value, -7.3, 1e-9);
}

TEST(TestTabGapPotential, TwoAndThreeBodyTermsGridMatching) {
    Species2Sorted pair_species("Fe", "Ni");
    ValuePtr<TwoBodyTransformation<2>> pair_trans = PairDistanceTransformation(CosCutoff(5.0, 1.0));
    SquaredExpKernel<1, 1> pair_kernel(1.0, std::array<Real, 1>{1.0});
    std::vector<Descriptor<2>> pair_sparse = {{2.0, 1.0}, {3.0, 1.0}};
    std::vector<Real> pair_coeffs = {0.5, -0.2};

    auto two_body_comp = TwoBodyGapComponent(pair_species, pair_trans, pair_kernel, pair_sparse, pair_coeffs);

    Species3AtomicSorted triplet_species("Fe", "Ni", "Fe");
    Angle3bTransformation triplet_trans(CosCutoff(4.0, 0.5));
    SquaredExpKernel<3, 1> triplet_kernel(1.0, std::array<Real, 3>{1.0, 1.0, 1.0});
    std::vector<Descriptor<4>> triplet_sparse = {{6.0, 0.0, 3.0, 1.0}};
    std::vector<Real> triplet_coeffs = {0.8};

    auto three_body_comp = ThreeBodyGapComponent<4, SquaredExpKernel<3, 1>>(
        triplet_species, ValuePtr<ThreeBodyTransformation<4>>(triplet_trans), triplet_kernel, triplet_sparse,
        triplet_coeffs
    );

    GapPotential gap;
    gap.addComponent(two_body_comp);
    gap.addComponent(three_body_comp);

    TabulationData tables = gap.tabulate({
        .max_cutoffs = gap.getCutoffs(),
        .n_grid_2b = 5000,
        .n_grid_3b = {20, 20, 20},
    });

    TabGapPotential tabgap(tables);

    // Verify 2-body spline predictions match table grid cell values
    const auto& two_body_map = tabgap.getTwoBodyComponents();
    ASSERT_TRUE(two_body_map.contains(pair_species));
    const auto& tg_2b = two_body_map.at(pair_species);
    const auto& grid_2b = tables.two_body_grids.value_grids.at(pair_species);

    for (const auto cell: grid_2b) {
        Real r = cell.pos[0];
        Real expected_val = cell.value;
        Real actual_val = tg_2b.getSpline()->interpolate(std::array{r}).value;
        EXPECT_NEAR(actual_val, expected_val, 1e-8);
    }

    // Verify 3-body spline predictions match table grid cell values
    const auto& three_body_map = tabgap.getThreeBodyComponents();
    ASSERT_TRUE(three_body_map.contains(triplet_species));
    const auto& tg_3b = three_body_map.at(triplet_species);
    const auto& grid_3b = tables.three_body_grids.value_grids.at(triplet_species);

    for (const auto cell: grid_3b) {
        Real expected_val = cell.value;
        Real actual_val = tg_3b.getSpline().interpolate(cell.pos).value;
        EXPECT_NEAR(actual_val, expected_val, 1e-8);
    }
}

TEST(TestTabGapPotential, EamPotentialTabulationAccuracy) {
    Species central_sp("Fe");
    PolycutoffPairFunction trans(4.0, 1.0, 2.0);

    auto eam_aggregator = TwoBodySum<1>("Fe");
    eam_aggregator.extend({"Fe", "Ni"}, trans);

    SquaredExpKernel<1, 0> kernel(1.0, std::array<Real, 1>{1.0});
    std::vector<Descriptor<1>> sparse_points = {{1.5}};
    std::vector<Real> coeffs = {1.2};

    auto eam_comp = ManyBodyGapComponent<1, SquaredExpKernel<1, 0>>(eam_aggregator, kernel, sparse_points, coeffs);

    GapPotential orig_gap;
    orig_gap.addComponent(eam_comp);

    TabulationData tables = orig_gap.tabulate({
        .max_cutoffs = orig_gap.getCutoffs(),
        .n_grid_2b = 500,
    });

    TabGapPotential tabgap(tables);

    Grid<3> grid_3b(
        std::array<size_t, 3>{10, 10, 10}, std::array<Real, 3>{0.3, 0.3, 0.2}, std::array<Real, 3>{0.5, 0.5, -1.0}
    );

    for (const auto cell: grid_3b) {
        Real r01 = cell.pos[0];
        Real r02 = cell.pos[1];
        Real cos12 = std::clamp(cell.pos[2], -1.0_r, 1.0_r);
        Real term = r01 * r01 + r02 * r02 - 2 * r01 * r02 * cos12;
        Real r12 = std::sqrt(std::max(0.0_r, term));

        Real sin12 = std::sqrt(std::max(0.0_r, 1.0_r - cos12 * cos12));

        Atoms triplet_atoms(
            {{0.0, 0.0, 0.0}, {r01, 0.0, 0.0}, {r02 * cos12, r02 * sin12, 0.0}},
            {central_sp, Species("Ni"), Species("Ni")}
        );

        auto orig_res = orig_gap.calculateEnergy(triplet_atoms);
        auto tabgap_res = tabgap.calculateEnergy(triplet_atoms);

        EXPECT_NEAR(tabgap_res.value, orig_res.value, 1e-6);

        ASSERT_EQ(tabgap_res.forces.size(), orig_res.forces.size());
        for (size_t i = 0; i < orig_res.forces.size(); ++i) {
            EXPECT_NEAR(tabgap_res.forces[i].x, orig_res.forces[i].x, 1e-3);
            EXPECT_NEAR(tabgap_res.forces[i].y, orig_res.forces[i].y, 1e-3);
            EXPECT_NEAR(tabgap_res.forces[i].z, orig_res.forces[i].z, 1e-3);
        }

        EXPECT_NEAR(tabgap_res.virials.xx, orig_res.virials.xx, 1e-3);
        EXPECT_NEAR(tabgap_res.virials.xy, orig_res.virials.xy, 1e-3);
        EXPECT_NEAR(tabgap_res.virials.xz, orig_res.virials.xz, 1e-3);
        EXPECT_NEAR(tabgap_res.virials.yz, orig_res.virials.yz, 1e-3);
        EXPECT_NEAR(tabgap_res.virials.yy, orig_res.virials.yy, 1e-3);
        EXPECT_NEAR(tabgap_res.virials.zz, orig_res.virials.zz, 1e-3);
    }
}

TEST(TestTabGapPotential, ThreeBodySameSpeciesClusterPermutationModeAndEnergyMatching) {
    Species3AtomicSorted same_species("Fe", "Fe", "Fe");
    Angle3bTransformation triplet_trans(CosCutoff(4.0, 0.5));
    SquaredExpKernel<3, 1> triplet_kernel(1.0, std::array<Real, 3>{1.0, 1.0, 1.0});
    std::vector<Descriptor<4>> triplet_sparse = {{5.0, 0.0, 2.5, 1.0}, {4.0, 0.5, 2.0, 1.0}};
    std::vector<Real> triplet_coeffs = {0.7, -0.3};

    auto three_body_comp = ThreeBodyGapComponent<4, SquaredExpKernel<3, 1>>(
        same_species, ValuePtr<ThreeBodyTransformation<4>>(triplet_trans), triplet_kernel, triplet_sparse,
        triplet_coeffs
    );

    GapPotential gap;
    gap.addComponent(three_body_comp);

    TabulationData tables = gap.tabulate({
        .max_cutoffs = gap.getCutoffs(),
        .n_grid_3b = {30, 30, 30},
    });

    TabGapPotential tabgap(tables);

    const auto& three_body_map = tabgap.getThreeBodyComponents();
    ASSERT_TRUE(three_body_map.contains(same_species));
    const auto& tg_3b = three_body_map.at(same_species);

    // Verify ClusterPermutationMode is NoNodePermutation
    EXPECT_EQ(tg_3b.getExpansion().getMode(), ClusterPermutationMode::NoNodePermutation);

    // Configuration of 4 Fe atoms (forming multiple triplets around central atoms)
    Atoms atoms(
        {
            {0.0, 0.0, 0.0},
            {2.0, 0.0, 0.0},
            {0.0, 2.2, 0.0},
            {1.5, 1.5, 0.0},
        },
        {Species("Fe"), Species("Fe"), Species("Fe"), Species("Fe")}
    );

    auto gap_res = gap.calculateEnergy(atoms);
    auto tabgap_res = tabgap.calculateEnergy(atoms);

    // If permutation mode was wrong (PermuteSameSpeciesNodes), tabgap_res would be 2x gap_res.
    EXPECT_NEAR(tabgap_res.value, gap_res.value, 1e-3);
    for (size_t i = 0; i < atoms.nAtoms(); ++i) {
        EXPECT_NEAR(tabgap_res.forces[i].x, gap_res.forces[i].x, 1e-2);
        EXPECT_NEAR(tabgap_res.forces[i].y, gap_res.forces[i].y, 1e-2);
        EXPECT_NEAR(tabgap_res.forces[i].z, gap_res.forces[i].z, 1e-2);
    }
}

TEST(TestTabGapPotential, ThreeBodySymmetricDistanceSwapCheck) {
    Species3AtomicSorted same_species("Fe", "Fe", "Fe");

    // Create an asymmetric 3D grid
    Grid<3> asymm_grid(
        std::array<size_t, 3>{4, 4, 4},
        std::array<Real, 3>{0.5, 0.5, 0.5},
        std::array<Real, 3>{1.0, 1.0, -1.0}
    );
    // Fill with values such that grid({1, 2, 0}) != grid({2, 1, 0})
    for (size_t i = 0; i < 4; ++i) {
        for (size_t j = 0; j < 4; ++j) {
            for (size_t k = 0; k < 4; ++k) {
                asymm_grid({i, j, k}) = static_cast<Real>(i * 10 + j + k);
            }
        }
    }

    CubicBSpline3D asymm_spline(asymm_grid);
    EXPECT_THROW(
        ThreeBodyTGComponent(same_species, asymm_spline),
        std::runtime_error
    );

    // Create a symmetric 3D grid: grid(i, j, k) == grid(j, i, k)
    Grid<3> symm_grid(
        std::array<size_t, 3>{4, 4, 4},
        std::array<Real, 3>{0.5, 0.5, 0.5},
        std::array<Real, 3>{1.0, 1.0, -1.0}
    );
    for (size_t i = 0; i < 4; ++i) {
        for (size_t j = 0; j < 4; ++j) {
            for (size_t k = 0; k < 4; ++k) {
                symm_grid({i, j, k}) = static_cast<Real>((i + j) * 10 + (i * j) + k);
            }
        }
    }

    CubicBSpline3D symm_spline(symm_grid);
    EXPECT_NO_THROW(
        ThreeBodyTGComponent(same_species, symm_spline)
    );
}
