#include <gtest/gtest.h>

#include "core/descriptors/TwoBodyDescriptor.hpp"

#include "core/neighbours/NeighbourFinder.hpp"
#include "utils/Utils.hpp"

using namespace jgap;

TEST(TestTwoBodyDescriptor, TrueUniformSimple) {

    auto structs = readXyz("test/resources/xyz-samples/iter-3-3-test.xyz");
    vector<AtomicStructure> selectedStructs;
    copy_n(structs.begin(), 20, back_inserter(selectedStructs));

    NeighbourFinder::findNeighbours(selectedStructs, 5.0);

    const auto params2b = TwoBodyDescriptorParams{
        .cutoff = 5.0,
        .kernelType = TwoBodyDescriptorParams::KernelType::GAUSS,
        .sparsificationMethod = TwoBodyDescriptorParams::SparsificationMethod::FULL_GRID_UNIFORM,
        .nSparsePointsPerSpeciesPair = 3,
        .energyScale = 10.0,
        .lengthScale = 1.0,
        .sparseRange = {1, 2},
        .speciesPairs = vector<SpeciesPair>{
            {"Mn", "Cr"},
            {"Cr", "Cr"}
        }
    };

    auto desc2b = TwoBodyDescriptor(params2b);
    desc2b.setSparsePoints(selectedStructs);
    const auto res = desc2b.serialize()["data"];
    ASSERT_EQ(res.size(), 2);
    ASSERT_EQ(res["Cr,Cr"]["sparse_points"].dump(), "[1.0,1.5,2.0]");
    ASSERT_EQ(res["Cr,Mn"]["sparse_points"].dump(), "[1.0,1.5,2.0]");
}

TEST(TestTwoBodyDescriptor, TrueUniformRangeFromFile) {

    auto structs = readXyz("test/resources/xyz-samples/iter-3-3-test.xyz");
    vector<AtomicStructure> selectedStructs;
    copy_n(structs.begin(), 20, back_inserter(selectedStructs));

    NeighbourFinder::findNeighbours(selectedStructs, 5.0);

    const auto params2b = TwoBodyDescriptorParams{
        .cutoff = 5.0,
        .kernelType = TwoBodyDescriptorParams::KernelType::GAUSS,
        .sparsificationMethod = TwoBodyDescriptorParams::SparsificationMethod::FULL_GRID_UNIFORM,
        .nSparsePointsPerSpeciesPair = 3,
        .energyScale = 10.0,
        .lengthScale = 1.0,
        .speciesPairs = vector<SpeciesPair>{
                {"Cr", "Cr"}
        }
    };

    auto desc2b = TwoBodyDescriptor(params2b);
    desc2b.setSparsePoints(selectedStructs);
    const auto res = desc2b.serialize()["data"];
    ASSERT_EQ(res.size(), 1);
    ASSERT_EQ(res["Cr,Cr"]["sparse_points"].dump(), "[3.48275862,3.68965517,3.89655172]");
}

TEST(TestTwoBodyDescriptor, TrueUniformSpeciesFromFile) {

    auto structs = readXyz("test/resources/xyz-samples/iter-3-3-test.xyz");
    vector<AtomicStructure> selectedStructs = structs;
    // copy_n(structs.begin(), 20, back_inserter(selectedStructs));

    NeighbourFinder::findNeighbours(selectedStructs, 5.0);

    const auto params2b = TwoBodyDescriptorParams{
        .cutoff = 5.0,
        .kernelType = TwoBodyDescriptorParams::KernelType::GAUSS,
        .sparsificationMethod = TwoBodyDescriptorParams::SparsificationMethod::FULL_GRID_UNIFORM,
        .nSparsePointsPerSpeciesPair = 3,
        .energyScale = 10.0,
        .lengthScale = 1.0
    };

    auto desc2b = TwoBodyDescriptor(params2b);
    desc2b.setSparsePoints(selectedStructs);
    cout << desc2b.serialize().dump()<<endl;
    const auto res = desc2b.serialize()["data"]; // num pairs
    ASSERT_EQ(res.size(), 4*5/2);
}

TEST(TestTwoBodyDescriptor, covariance) {

    auto structs = readXyz("test/resources/xyz-samples/iter-3-3-test.xyz");
    NeighbourFinder::findNeighbours(structs, 5.0);
    vector<AtomicStructure> selectedStructs = structs;
    // copy_n(structs.begin(), 20, back_inserter(selectedStructs));

    const auto params2b = TwoBodyDescriptorParams{
        .cutoff = 3.5, // => no self-interaction in x direction
        .cutoffTransitionWidth = 0.7, // => should alter CrCr 2.89132 -> ~9.94112 to ~9.52948
        .kernelType = TwoBodyDescriptorParams::KernelType::GAUSS,
        .sparsificationMethod = TwoBodyDescriptorParams::SparsificationMethod::FULL_GRID_UNIFORM,
        .nSparsePointsPerSpeciesPair = 2,
        .energyScale = 10.0,
        .lengthScale = 1.0,
        .sparseRange = {2, 3},
        .speciesPairs = vector<SpeciesPair>{
            {"Mn", "Cr"},
            {"Cr", "Cr"}
        }
    };

    auto desc2b = TwoBodyDescriptor(params2b);
    desc2b.setSparsePoints(selectedStructs);

    const auto resSparse = desc2b.serialize()["data"];
    ASSERT_EQ(resSparse.size(), 2);
    ASSERT_EQ(resSparse["Cr,Cr"]["sparse_points"].dump(), "[2.0,3.0]");
    ASSERT_EQ(resSparse["Cr,Mn"]["sparse_points"].dump(), "[2.0,3.0]");

    auto covariance = desc2b.covariate(
        structs[128] // 4 atoms: crmnfeni elongated in x(>4) / ~2.7 in z and y
    );
    ASSERT_EQ(covariance.size(), 4);

    auto U = vector{
        291, 378, // self interact
        725, 677 // do not self interact
    }; // => double chack changes if this fails
    for (int i = 0; i < 4; i++) EXPECT_NEAR(covariance[i].total, U[i], 2);

    // Mn atom
    // no CrCr covariance
    for (int j = 0; j < 2; j++) {
        EXPECT_NEAR(covariance[j].derivatives[0].x, 0, 1e-9);
        EXPECT_NEAR(covariance[j].derivatives[0].y, 0, 1e-9);
        EXPECT_NEAR(covariance[j].derivatives[0].z, 0, 1e-9);
    }
    // sparseDist = 2 & (x++ => dist->2) => K_se++
    ASSERT_TRUE(covariance[2].derivatives[0].x > 0);
    // opposite
    ASSERT_TRUE(covariance[3].derivatives[0].x < 0);

    // Cr atom
    // pull opposite to manganese
    ASSERT_TRUE(covariance[2].derivatives[1].x < 0);
    ASSERT_TRUE(covariance[3].derivatives[1].x > 0);
    // symmetry
    for (int j = 0; j < 2; j++) {
        EXPECT_NEAR(covariance[j].derivatives[1].x, 0.0, 1e-9);
        EXPECT_NEAR(covariance[j].derivatives[1].y, 0.0, 1e-9);
        EXPECT_NEAR(covariance[j].derivatives[1].z, 0.0, 1e-9);
    }

    // Non-Cr/Mn atoms
    for (int j = 0; j < 4; j++) {
        for (int i = 2; i < 3; i++) EXPECT_NEAR(covariance[j].derivatives[i].x, 0, 1e-9);
        for (int i = 2; i < 3; i++) EXPECT_NEAR(covariance[j].derivatives[i].y, 0, 1e-9);
        for (int i = 2; i < 3; i++) EXPECT_NEAR(covariance[j].derivatives[i].z, 0, 1e-9);
    }
}

AtomicStructure makeDimer(Species species, double x) {
    return AtomicStructure{
        .atoms = {
            AtomData{
                .species = std::move(species),
                .position = {0,0,0}
            },
            AtomData{
                .species = std::move(species),
                .position = {x,0,0}
            }
        },
        .latticeVectors = {
            Vector3{30.0, 0.0, 0.0},
            Vector3{0.0, 30.0, 0.0},
            Vector3{0.0, 0.0, 30.0},
        }
    };
}

TEST(TestTwoBodyDescriptor, dimerRepSign) {
    vector structures = {
        makeDimer("Fe", 1.0),
        makeDimer("Fe", 2.5),
        makeDimer("Fe", 4),
        };
    NeighbourFinder::findNeighbours(structures, 5.0);

    const auto params2b = TwoBodyDescriptorParams{
        .cutoff = 5, // => no self-interaction in x direction
        .cutoffTransitionWidth = 0.7, // => should alter CrCr 2.89132 -> ~9.94112 to ~9.52948
        .kernelType = TwoBodyDescriptorParams::KernelType::GAUSS,
        .sparsificationMethod = TwoBodyDescriptorParams::SparsificationMethod::FULL_GRID_UNIFORM,
        .nSparsePointsPerSpeciesPair = 2,
        .energyScale = 1.0,
        .lengthScale = 1.0,
        .sparseRange = {2, 3}
    };

    auto desc2b = TwoBodyDescriptor(params2b);
    desc2b.setSparsePoints(structures);

    auto res0 = desc2b.covariate(structures[0]);
    ASSERT_TRUE(res0[0].derivatives[0].x < 0);
    ASSERT_TRUE(res0[1].derivatives[0].x < 0);
    ASSERT_TRUE(res0[0].derivatives[1].x > 0);
    ASSERT_TRUE(res0[1].derivatives[1].x > 0);

    auto res1 = desc2b.covariate(structures[1]);
    ASSERT_TRUE(res1[0].derivatives[0].x > 0);
    ASSERT_TRUE(res1[1].derivatives[0].x < 0);
    ASSERT_TRUE(res1[0].derivatives[1].x < 0);
    ASSERT_TRUE(res1[1].derivatives[1].x > 0);

    auto res2 = desc2b.covariate(structures[2]);
    ASSERT_TRUE(res2[0].derivatives[0].x > 0);
    ASSERT_TRUE(res2[1].derivatives[0].x > 0);
    ASSERT_TRUE(res2[0].derivatives[1].x < 0);
    ASSERT_TRUE(res2[1].derivatives[1].x < 0);
}


auto equilateralTriangle_2b = AtomicStructure{
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

TwoBodyDescriptor setupDesc2bSimple_2b() {
    const auto params2b = TwoBodyDescriptorParams{
        .cutoff = 10, // => no self-interaction in x direction
        .cutoffTransitionWidth = 0.7, // => should alter CrCr 2.89132 -> ~9.94112 to ~9.52948
        .kernelType = TwoBodyDescriptorParams::KernelType::GAUSS,
        .sparsificationMethod = TwoBodyDescriptorParams::SparsificationMethod::FULL_GRID_UNIFORM,
        .nSparsePointsPerSpeciesPair = 2,
        .energyScale = 1.0,
        .lengthScale = 1.0,
        .sparseRange = {2, 3}
    };

    auto desc2b = TwoBodyDescriptor(params2b);
    desc2b.setSparsePoints({
        {SpeciesPair{"Fe", "Fe"}, vector{2.0, 4.0}}
    });

    return desc2b;
}

TEST(TestTwoBodyDescriptor, twoBodyEquilateralTriangle) {
    auto desc = setupDesc2bSimple_2b();
    NeighbourFinder::findNeighbours(equilateralTriangle_2b, 5.0);
    auto res = desc.covariate(equilateralTriangle_2b);

    // energy
    ASSERT_NEAR(res[0].total, 2*3*exp(-0.5), 2e-4);
    ASSERT_NEAR(res[1].total, 2*3*exp(-0.5), 2e-4);

    // derivatives
    double dkdr_at_2 = -2.0 * exp(-0.5);
    double dkdr_at_4 = +2.0 * exp(-0.5);

    auto atoms = equilateralTriangle_2b.atoms;

    for (int i: {0, 1, 2}) {
        Vector3 derivatives_at2{0,0,0}, derivatives_at4{0,0,0};
        for (int j: {0, 1, 2}) {
            if (i == j) continue;
            derivatives_at2 = derivatives_at2 + (atoms[i].position - atoms[j].position).normalize() * dkdr_at_2;
            derivatives_at4 = derivatives_at4 + (atoms[i].position - atoms[j].position).normalize() * dkdr_at_4;
        }
        cout << i << " at 2: " << derivatives_at2.toString() << endl;
        cout << i << " at 4: " << derivatives_at4.toString() << endl;
        ASSERT_NEAR((res[0].derivatives[i]-derivatives_at2).norm(), 0, 1e-4);
        ASSERT_NEAR((res[1].derivatives[i]-derivatives_at4).norm(), 0, 1e-4);
    }

    desc.serialize();
}

TEST(TestTwoBodyDescriptor, doubleBoxDoubleEnergy) {
    vector structs = readXyz("test/resources/xyz-samples/FeOnly.xyz");
    NeighbourFinder::findNeighbours(structs, 5.0);

    const auto params2b = TwoBodyDescriptorParams{
        .cutoff = 5.0,
        .kernelType = TwoBodyDescriptorParams::KernelType::GAUSS,
        .sparsificationMethod = TwoBodyDescriptorParams::SparsificationMethod::FULL_GRID_UNIFORM,
        .nSparsePointsPerSpeciesPair = 3,
        .energyScale = 10.0,
        .lengthScale = 1.0,
        .sparseRange = {1, 2}
    };

    auto desc2b = TwoBodyDescriptor(params2b);
    desc2b.setSparsePoints(structs);
    desc2b.setCoefficients({1, 1, 1});

    for (size_t i = 71; i < 440; i++) {
        auto structure = structs[i];
        cout << structure.atoms.size() << endl;
        auto rep = structure.repeat(2,2,2);
        NeighbourFinder::findNeighbours(rep, 5.0);
        auto predOrig = desc2b.predict(structure);
        auto predRep = desc2b.predict(rep);
        ASSERT_NEAR(
            predOrig.energy.value() * 8.0,
            predRep.energy.value(),
            1e-4
            );
        Vector3 forceSumOrig{0,0,0};
        for (auto f: predOrig.forces.value()) {
            forceSumOrig = forceSumOrig + f;
        }
        Vector3 forceSumRep{0,0,0};
        for (auto f: predRep.forces.value()) {
            forceSumRep = forceSumRep + f;
        }
        ASSERT_NEAR(forceSumOrig.norm() * 8, forceSumRep.norm(), 1e-4);
    }
}