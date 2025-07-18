#include <gtest/gtest.h>

#include "core/descriptors/ThreeBodyDescriptor.hpp"

#include "core/neighbours/NeighbourFinder.hpp"
#include "utils/Utils.hpp"

using namespace jgap;

AtomicStructure pythagorian3b;
nlohmann::json params3b;

void initPythagorian3b() {
    pythagorian3b = AtomicStructure{
        .lattice = {
            Vector3{20, 0, 0},
            Vector3{0, 20, 0},
            Vector3{0, 0, 20},
        },
        .atoms = {
            AtomData{
                .position = Vector3{0, 0, 0},
                .species = "Fe"
            },
            AtomData{
                .position = Vector3{4, 0, 0},
                .species = "Fe"
            },
            AtomData{
                .position = Vector3{0, 3, 0},
                .species = "Fe"
            }
        }
    };
    NeighbourFinder::findNeighbours(pythagorian3b, 10.0);
}

void initThreeBodyParams() {
    params3b = nlohmann::json::parse(R"(
    {
        "kernel": {
            "type": "squared_exp",
            "length_scale": 1.0,
            "energy_scale": 1.0,
            "cutoff": {
                "type": "perriot",
                "cutoff_transition_width": 0.5,
                "cutoff": 10.0
            }
        },
        "sparse_data": {}
    }
    )");
}

TEST(TestThreeBodyDescriptor, CovarianceEnergy) {
    initPythagorian3b();
    initThreeBodyParams();

    auto desc3b = ThreeBodyDescriptor(params3b);
    desc3b.setSparsePoints(map<SpeciesTriplet, vector<Vector3>>{
        {
            SpeciesTriplet{.root = "Fe", .nodes = SpeciesPair{"Fe", "Fe"}}, {Vector3{8, 2, 4}}
        }
    });

    auto result = desc3b.covariate({pythagorian3b});
    ASSERT_EQ(result.size(), 1);
    ASSERT_NEAR(result[0].total, 2.0 * (2.0 * exp(-3.0/2.0) + exp(-2.0)), 1e-6);
}

TEST(TestThreeBodyDescriptor, CovarianceDerivatives) {
    initPythagorian3b();
    pythagorian3b.atoms[1].species = "Ni";
    pythagorian3b.atoms[2].species = "Cr";

    initThreeBodyParams();

    auto desc3b = ThreeBodyDescriptor(params3b);
    desc3b.setSparsePoints(map<SpeciesTriplet, vector<Vector3>>{
        {
            SpeciesTriplet{.root = "Fe", .nodes = SpeciesPair{"Ni", "Cr"}},
            {Vector3{8, 2, 4}}
        }
    });

    auto result = desc3b.covariate({pythagorian3b});
    ASSERT_EQ(result.size(), 1);
    ASSERT_NEAR(result[0].total, 2.0 * exp(-3.0/2.0), 1e-6);
    ASSERT_NEAR(result[0].forces[0].z, 0, 1e-6);
    ASSERT_NEAR(result[0].forces[1].z, 0, 1e-6);
    ASSERT_NEAR(result[0].forces[2].z, 0, 1e-6);

    ASSERT_NEAR(result[0].forces[0].x, 2.0 * -3.0 * exp(-3.0/2.0), 1e-6);
    ASSERT_NEAR(result[0].forces[0].y, 2.0 * 1.0 * exp(-3.0/2.0), 1e-6);

    ASSERT_NEAR(result[0].forces[1].x, 2.0 * 2.2 * exp(-3.0/2.0), 1e-6);
    ASSERT_NEAR(result[0].forces[1].y, 2.0 * 0.6 * exp(-3.0/2.0), 1e-6);

    ASSERT_NEAR(result[0].forces[2].x, 2.0 * 0.8 * exp(-3.0/2.0), 1e-6);
    ASSERT_NEAR(result[0].forces[2].y, 2.0 * -1.6 * exp(-3.0/2.0), 1e-6);
}

