#include <algorithm>
#include <filesystem>
#include <numeric>
#include <random>
#include <gtest/gtest.h>

#include "jgap/core/atomic/Atoms.hpp"
#include "jgap/core/potentials/gap/GapPotential.hpp"
#include "jgap/serialization/SerializationRegistry.hpp"

using namespace jgap;

TEST(TestGapPotential, EmptyPotential) {
    GapPotential pot;
    EXPECT_EQ(pot.getComponents().size(), 0);
    EXPECT_DOUBLE_EQ(pot.getCutoffs().maxOverall(), 0.0);

    Atoms atoms({{0, 0, 0}, {1, 1, 1}}, {Species("Fe"), Species("Fe")});
    auto res = pot.calculateEnergy(atoms);
    EXPECT_DOUBLE_EQ(res.value, 0.0);
    EXPECT_EQ(res.forces.size(), 2);
    EXPECT_DOUBLE_EQ(res.forces[0].x, 0.0);
}

TEST(TestGapPotential, QuipReferencePredictionSubset) {
    namespace fs = std::filesystem;
    fs::path h5_path = "test/resources/reference/reference_pots/feni_200_100_0/gap.h5";
    fs::path pred_xyz = "test/resources/reference/reference_preds/feni_200_100_0/pred.xyz";
    ASSERT_TRUE(fs::exists(h5_path)) << "Reference data not found: " << h5_path;
    ASSERT_TRUE(fs::exists(pred_xyz)) << "Reference data not found: " << pred_xyz;

    ValuePtr<Potential> potential = SerializationRegistry<Potential>::deserialize(h5_path.string());
    ASSERT_NE(potential.get(), nullptr);

    MainXYZPropertyNames prop_names;
    prop_names.virials = "virial";
    auto ref_frames = Atoms::readAtoms(pred_xyz.string(), prop_names);
    ASSERT_FALSE(ref_frames.empty());

    std::vector<size_t> indices(ref_frames.size());
    std::iota(indices.begin(), indices.end(), 0);
    std::mt19937 rng(42);
    std::shuffle(indices.begin(), indices.end(), rng);
    size_t n_to_check = std::min<size_t>(100, indices.size());

    constexpr double meV = 1e3;
    for (size_t i = 0; i < n_to_check; ++i) {
        size_t f = indices[i];
        const auto& ref_atoms = ref_frames[f];
        size_t n = ref_atoms.nAtoms();
        ASSERT_GT(n, 0);

        Atoms test_atoms = ref_atoms;
        auto res = potential->calculateEnergy(test_atoms);

        double e_jgap = res.value / static_cast<double>(n) * meV;
        double e_ref = ref_atoms.getEnergy().value_or(0.0) / static_cast<double>(n) * meV;
        EXPECT_NEAR(e_jgap, e_ref, 0.05);

        const auto ref_forces = ref_atoms.getForces();
        if (ref_forces.has_value()) {
            ASSERT_EQ(res.forces.size(), ref_forces->size());
            for (size_t a = 0; a < res.forces.size(); ++a) {
                EXPECT_NEAR(res.forces[a].x * meV, (*ref_forces)[a].x * meV, 0.5);
                EXPECT_NEAR(res.forces[a].y * meV, (*ref_forces)[a].y * meV, 0.5);
                EXPECT_NEAR(res.forces[a].z * meV, (*ref_forces)[a].z * meV, 0.5);
            }
        }
    }
}
