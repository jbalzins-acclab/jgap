#include <gtest/gtest.h>
#include <gtest/gtest.h>

#include "core/descriptors/EamDescriptor.hpp"
#include "core/descriptors/eam/EamDensityCalculator.hpp"
#include "core/descriptors/eam/pair_functions/PolycutoffPairFunction.hpp"
#include "core/neighbours/NeighbourFinder.hpp"
#include "data/BasicDataTypes.hpp"

using namespace jgap;

auto equilateralTriangleEAM = AtomicStructure{
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
            .position = {3.0, 0.0, 0.0},
            .species = "Fe"
        },
        AtomData{
            .position = {1.5, 2.598, 0.0},
            .species = "Fe"
        }
    }
};

TEST(EamDescriptorTest, equilateralTriangle) {
    auto params = EamDescriptorParams{
        .defaultDensityCalculationParams = EamDensityCalculationParams{
            .pairFunction = EamDensityCalculationParams::PairFunction::POLYCUTOFF,
            .speciesUse = EamDensityCalculationParams::SpeciesUse::BLIND,
            .cutoff = 6,
            .rmin = 0.0
            },
        .nSparsePoints = 2,
        .energyScale = 1,
        .lengthScale = 1,
        .sparseRange = {0, 2},
        .kernelType = EamDescriptorParams::KernelType::GAUSS,
        .sparsificationMethod = EamDescriptorParams::SparsificationMethod::FULL_GRID_UNIFORM
    };

    auto desc = EamDescriptor(params);

    NeighbourFinder::findNeighbours(equilateralTriangleEAM, 6.0);
    desc.setSparsePoints({equilateralTriangleEAM});

    auto result = desc.covariate(equilateralTriangleEAM);

    ASSERT_NEAR(result[0].total, 3*exp(-0.5), 1e-4);
    ASSERT_NEAR(result[1].total, 3*exp(-0.5), 1e-4);

    auto atoms = equilateralTriangleEAM.atoms;
    for (int i: {0, 1, 2}) {
        Vector3 total0{0,0,0}, total2{0,0,0};
        for (int j: {0, 1, 2}) {
            if (i == j) continue;
            total0 = total0 + (atoms[i].position - atoms[j].position).normalize() * 0.3125 * exp(-0.5);
            total2 = total2 + (atoms[i].position - atoms[j].position).normalize() * -0.3125 * exp(-0.5);
        }
        ASSERT_NEAR((result[0].derivatives[i] - total0).norm(), 0, 1e-4);
        ASSERT_NEAR((result[1].derivatives[i] - total2).norm(), 0, 1e-4);
    }
}