#include <gtest/gtest.h>

#include "core/descriptors/EamDescriptor.hpp"
#include "core/neighbours/NeighbourFinder.hpp"
#include "data/BasicDataTypes.hpp"

using namespace jgap;

TEST(EamDescriptorTest, equilateralTriangleEamTest) {

    auto equilateralTriangleEAM = AtomicStructure{
        .lattice = {
            Vector3{100.0, 0.0, 0.0},
            Vector3{0.0, 100.0, 0.0},
            Vector3{0.0, 0.0, 100.0}
        },
        .atoms = {
            AtomData{
                .position = {0.0, 0.0, 0.0},
                .species = "Fe"
            },
            AtomData{
                .position = {3.0, 0.0, 0.0},
                .species = "Fe"
            },
            AtomData{
                .position = {1.5, 2.598, 0.0},
                .species = "Fe"
            }
        }
    };

    auto params = nlohmann::json::parse(R"(
    {
        "kernel": {
            "type": "squared_exp",
            "length_scale": 1.0,
            "energy_scale": 1.0
        },
        "pair_functions": [
            {
                "type": "polycutoff",
                "r_min": 0,
                "cutoff": 6.0
            }
        ],
        "sparse_data": {
            "Fe": {
                "sparse_points": [0, 2]
            }
        }
    }
    )");
    auto desc = EamDescriptor(params);

    NeighbourFinder::findNeighbours(equilateralTriangleEAM, 6.0);
    // desc.setSparsePoints({equilateralTriangleEAM});

    auto result = desc.covariate(equilateralTriangleEAM);

    ASSERT_NEAR(result[0].total, 3.0*exp(-0.5), 1e-4);
    ASSERT_NEAR(result[1].total, 3.0*exp(-0.5), 1e-4);

    auto atoms = equilateralTriangleEAM.atoms;
    for (int i: {0, 1, 2}) {
        Vector3 total0{0,0,0}, total2{0,0,0};
        for (int j: {0, 1, 2}) {
            if (i == j) continue;
            total0 = total0 + (atoms[i].position - atoms[j].position).normalize() * 0.3125 * exp(-0.5) * 2.0;
            total2 = total2 + (atoms[i].position - atoms[j].position).normalize() * -0.3125 * exp(-0.5) * 2.0;
        }
        ASSERT_NEAR((result[0].derivatives[i] - total0).norm(), 0, 1e-4);
        ASSERT_NEAR((result[1].derivatives[i] - total2).norm(), 0, 1e-4);
    }
}