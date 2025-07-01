#include <gtest/gtest.h>

#include "core/descriptors/eam/EamDensityCalculator.hpp"
#include "core/descriptors/eam/pair_functions/PolycutoffPairFunction.hpp"
#include "core/neighbours/NeighbourFinder.hpp"
#include "data/BasicDataTypes.hpp"

using namespace jgap;

auto equilateralTriangleEDC = AtomicStructure{
    .latticeVectors = {
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
            .position = {3, 0.0, 0.0},
            .species = "Fe"
        },
        AtomData{
            .position = {1.5, 2.598, 0.0},
            .species = "Fe"
        }
    }
};

TEST(TestEamDensityCalc, equilateralTriangle) {

    auto params = EamDescriptorParams{
        .defaultDensityCalculationParams = EamDensityCalculationParams{
            .pairFunction = EamDensityCalculationParams::PairFunction::POLYCUTOFF,
            .speciesUse = EamDensityCalculationParams::SpeciesUse::BLIND,
            .cutoff = 6,
            .rmin = 0.0
            }
    };
    EamDensityCalculator calc(params);
    NeighbourFinder::findNeighbours(equilateralTriangleEDC, 6.0);

    auto res = calc.calculate(equilateralTriangleEDC);
    ASSERT_EQ(res.size(), 3);
    for (int i: {0, 1, 2}) {
        ASSERT_NEAR(res[i].density, 1.0, 1e-4);
        ASSERT_EQ(res[0].densityDerivatives.size(), 2);
        for (int j: {0, 1}) {
            ASSERT_NEAR(res[0].densityDerivatives[j].second, -0.3125, 1e-4);
        }
    }
}