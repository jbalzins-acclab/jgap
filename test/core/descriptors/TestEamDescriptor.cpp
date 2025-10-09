#include <gtest/gtest.h>

#include "core/descriptors/eam/EamDescriptor.hpp"
#include "core/neighbours/NeighbourFinder.hpp"
#include "data/BasicDataTypes.hpp"
#include "ParserRegistryAuto.hpp"

using namespace jgap;

TEST(EamDescriptorTest, equilateralTriangleEamTest) {

    const vector positions = {
        Vector3{0.0, 0.0, 0.0},
        Vector3{3.0, 0.0, 0.0},
        Vector3{1.5, 2.598, 0.0}
    };
    auto equilateralTriangleEAM = AtomicStructure{
        .lattice = {
            Vector3{100.0, 0.0, 0.0},
            Vector3{0.0, 100.0, 0.0},
            Vector3{0.0, 0.0, 100.0}
        },
        .positions = positions,
        .species = {"Fe", "Fe", "Fe"}
    };

    const auto params = nlohmann::json::parse(R"(
    {
        "kernels": [
            {
                "species": "Fe",
                "type": "squared_exp",
                "length_scale": 1.0,
                "energy_scale": 1.0,
                "density": 0
            },
            {
                "species": "Fe",
                "type": "squared_exp",
                "length_scale": 1.0,
                "energy_scale": 1.0,
                "density": 2
            }
        ],
        "pair_functions": [
            {
                "type": "polycutoff",
                "r_min": 0,
                "cutoff": 6.0
            }
        ]
    }
    )");
    auto desc = EamDescriptor(params);

    NeighbourFinder::findNeighbours(equilateralTriangleEAM, 6.0);
    // desc.setSparsePoints({equilateralTriangleEAM});

    auto result = desc.covariate(equilateralTriangleEAM);

    /*
    for sanity:

    sparse point 0:
    cout << result[0].forces[0].toString() << endl;
    cout << result[0].forces[1].toString() << endl;
    cout << result[0].forces[2].toString() << endl;
    OUT:
    0.568627, 0.328292, 0.000000
    -0.568627, 0.328292, 0.000000
    0.000000, -0.656584, 0.000000

    -> decrease separation => increase density => covariance K(\rho, 0) decreases
    -> F should point toward direction where K would decrease => F is attractive
    */
    ASSERT_NEAR(result[0].total, 3.0*exp(-0.5), 1e-4);
    ASSERT_NEAR(result[1].total, 3.0*exp(-0.5), 1e-4);

    for (int i: {0, 1, 2}) {
        Vector3 total0{0,0,0}, total2{0,0,0};
        for (int j: {0, 1, 2}) {
            if (i == j) continue;
            total0 = total0 + (positions[j] - positions[i]).normalize() * 0.3125 * exp(-0.5) * 2.0;
            total2 = total2 + (positions[j] - positions[i]).normalize() * -0.3125 * exp(-0.5) * 2.0;
        }
        ASSERT_NEAR((result[0].forces[i] - total0).len(), 0, 1e-4);
        ASSERT_NEAR((result[1].forces[i] - total2).len(), 0, 1e-4);
    }
}