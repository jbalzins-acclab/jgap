#include <gtest/gtest.h>

#include "core/descriptors/ThreeBodyDescriptor.hpp"

#include "core/neighbours/NeighbourFinder.hpp"
#include "data/params/TwoBodyDescriptorParams.hpp"
#include "utils/Utils.hpp"

using namespace jgap;

AtomicStructure pythagorian3b;
ThreeBodyDescriptorParams params3b;

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
    params3b = ThreeBodyDescriptorParams{
        .cutoff = 10.0,
        .kernelType = ThreeBodyDescriptorParams::KernelType::GAUSS,
        .sparsificationMethod = ThreeBodyDescriptorParams::SparsificationMethod::FULL_GRID_UNIFORM,
        .nSparsePointsPerSpeciesPerDirection = array<size_t, 3>{2,2,2},
        .energyScale = 1.0,
        .lengthScale = 1.0
    };
}

TEST(TestThreeBodyDescriptor, TrueUniformSimple) {
    initPythagorian3b();
    initThreeBodyParams();

    auto desc3b = ThreeBodyDescriptor(params3b);
    desc3b.setSparsePoints({pythagorian3b});
    // cout << desc3b.serialize();

    auto result = desc3b.serialize();

    ASSERT_EQ(result["data"].size(), 1);
    ASSERT_EQ(result["data"]["Fe,Fe,Fe"]["sparse_points"].size(), 8);

    ASSERT_NEAR(result["data"]["Fe,Fe,Fe"]["sparse_points"][0]["x"].get<double>(), 7.0, 1e-6);
    ASSERT_NEAR(result["data"]["Fe,Fe,Fe"]["sparse_points"][0]["y"].get<double>(), 1.0, 1e-6);
    ASSERT_NEAR(result["data"]["Fe,Fe,Fe"]["sparse_points"][0]["z"].get<double>(), 3.0, 1e-6);

    ASSERT_NEAR(result["data"]["Fe,Fe,Fe"]["sparse_points"][7]["x"].get<double>(), 9.0, 1e-6);
    ASSERT_NEAR(result["data"]["Fe,Fe,Fe"]["sparse_points"][7]["y"].get<double>(), 4.0, 1e-6);
    ASSERT_NEAR(result["data"]["Fe,Fe,Fe"]["sparse_points"][7]["z"].get<double>(), 5.0, 1e-6);
}

TEST(TestThreeBodyDescriptor, CovarianceEnergy) {
    initPythagorian3b();
    initThreeBodyParams();

    auto desc3b = ThreeBodyDescriptor(params3b);
    desc3b.setSparsePoints(map<SpeciesTriplet, vector<Vector3>>{
        {
            SpeciesTriplet{.root = "Fe", .nodes = SpeciesPair{"Fe", "Fe"}},
            {Vector3{8, 2, 4}}
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
    ASSERT_NEAR(result[0].derivatives[0].z, 0, 1e-6);
    ASSERT_NEAR(result[0].derivatives[1].z, 0, 1e-6);
    ASSERT_NEAR(result[0].derivatives[2].z, 0, 1e-6);

    ASSERT_NEAR(result[0].derivatives[0].x, 2.0 * -3.0 * exp(-3.0/2.0), 1e-6);
    ASSERT_NEAR(result[0].derivatives[0].y, 2.0 * 1.0 * exp(-3.0/2.0), 1e-6);

    ASSERT_NEAR(result[0].derivatives[1].x, 2.0 * 2.2 * exp(-3.0/2.0), 1e-6);
    ASSERT_NEAR(result[0].derivatives[1].y, 2.0 * 0.6 * exp(-3.0/2.0), 1e-6);

    ASSERT_NEAR(result[0].derivatives[2].x, 2.0 * 0.8 * exp(-3.0/2.0), 1e-6);
    ASSERT_NEAR(result[0].derivatives[2].y, 2.0 * -1.6 * exp(-3.0/2.0), 1e-6);
}

