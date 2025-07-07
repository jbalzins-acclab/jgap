#include <gtest/gtest.h>

#include "core/InRamJgapFit.hpp"
#include "core/descriptors/EamDescriptor.hpp"
#include "core/descriptors/ThreeBodyDescriptor.hpp"
#include "core/descriptors/TwoBodyDescriptor.hpp"
#include "core/matrices/sigmas/SimpleSigmaRules.hpp"
#include "core/neighbours/NeighbourFinder.hpp"
#include "data/BasicDataTypes.hpp"
#include "data/params/EamDescriptorParams.hpp"
#include "data/params/ThreeBodyDescriptorParams.hpp"
#include "utils/Utils.hpp"

using namespace jgap;

AtomicStructure equilateralTriangle;
void setupEquilateralTriangle() {
    equilateralTriangle = AtomicStructure{
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
}

TEST(TestInRamJgap, twoBodyEquilateralTriangleAtEquilibriumQuipCompatibility) {
    setupEquilateralTriangle();
    const auto params2b = TwoBodyDescriptorParams{
        .cutoff = 10,
        .cutoffTransitionWidth = 0.7,
        .kernelType = TwoBodyDescriptorParams::KernelType::GAUSS,
        .sparsificationMethod = TwoBodyDescriptorParams::SparsificationMethod::FULL_GRID_UNIFORM,
        .nSparsePointsPerSpeciesPair = 1,
        .energyScale = 1.0,
        .lengthScale = 1.0,
        .sparseRange = {2, 3}
    };

    auto desc2b = TwoBodyDescriptor(params2b);
    desc2b.setSparsePoints({
        {SpeciesPair{"Fe", "Fe"}, vector{3.0}}
    });

    auto fit = InRamJgapFit(
        {make_shared<TwoBodyDescriptor>(desc2b)},
        {},
        make_shared<SimpleSigmaRules>(3, 10, 5, 5),
        1e-8
        );

    equilateralTriangle.energy = 1;
    equilateralTriangle.atoms[0].force = Vector3{0.0, 0.0, 0.0};
    equilateralTriangle.atoms[1].force = Vector3{0.0, 0.0, 0.0};
    equilateralTriangle.atoms[2].force = Vector3{0.0, 0.0, 0.0};
    //NeighbourFinder::findNeighbours(equilateralTriangle, 10.0);
    //auto a = desc2b.covariate(equilateralTriangle);
    auto pot = fit.fitGAP(vector{equilateralTriangle});

    auto coeffs = pot->serialize()["descriptors"][0]["data"]["Fe,Fe"]["coefficients"];
    cout << coeffs.dump() << endl;
    auto c = coeffs[0].get<double>();
    /*
     * if a - total covariance
     * => a = 2(idk why - ?counted from both ends?) * 3 * (K(3, 3) = 1) = 6
     * => c = a / (\sigma^2_{total-i.e. per all atoms} + a^2)
     * => variance + a^2 = a/c
     * => variance = a/c - a^2
     * \sigma / "per atom" = 3 => total variance = 3^2 * 3(atoms)
     */
    auto totalVariance = 6.0 / c - 36.0;
    ASSERT_NEAR(totalVariance, 27, 1e-6);
    // GPWRITE=1 VERBOSE=1 gap_fit at_file=train.xyz e0=0.0 sparse_jitter=1e-8 do_copy_at_file=False gp_file=gap_out.xml rnd_seed=999 default_sigma="{4.0 1.0 0 0}" gap="{distance_2b cutoff=10.0 covariance_type=ard_se delta=1.0 theta_uniform=1.0 sparse_method=FILE sparse_file=mysparse.desc print_sparse_index=sparse_indices_2b.out}"
}

TEST(TestInRamJgap, twoAtomsWithForceQuipCompatibility) {
    // TODO: split
    auto twoAtoms = AtomicStructure{
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
            }
        }
    };

    const auto params2b = TwoBodyDescriptorParams{
        .cutoff = 10,
        .cutoffTransitionWidth = 0.7,
        .kernelType = TwoBodyDescriptorParams::KernelType::GAUSS,
        .sparsificationMethod = TwoBodyDescriptorParams::SparsificationMethod::FULL_GRID_UNIFORM,
        .nSparsePointsPerSpeciesPair = 1,
        .energyScale = 1.0,
        .lengthScale = 1.0,
        .sparseRange = {2, 3}
    };

    auto desc2b = TwoBodyDescriptor(params2b);
    desc2b.setSparsePoints({
        {SpeciesPair{"Fe", "Fe"}, vector{2.0,4.0}}
    });

    auto fit = InRamJgapFit(
        {make_shared<TwoBodyDescriptor>(desc2b)},
        {},
        make_shared<SimpleSigmaRules>(1, 1, 5, 5),
        1e-8
        );

    twoAtoms.energy = 1;
    twoAtoms.atoms[0].force = Vector3{0.0, 0.0, 0.0};
    twoAtoms.atoms[1].force = Vector3{0.0, 0.0, 0.0};
    auto pot = fit.fitGAP(vector{twoAtoms});

    auto coeffs = pot->serialize()["descriptors"][0]["data"]["Fe,Fe"]["coefficients"];
    cout << coeffs.dump() << endl;
    auto a0 = coeffs[0].get<double>(), a1 = coeffs[1].get<double>();
    ASSERT_NEAR(a0, 0.23266775926220731, 1e-6);
    ASSERT_NEAR(a1, 0.23266775926220731, 1e-6);

    //////////////////////////
    auto desc2b2 = TwoBodyDescriptor(params2b);
    desc2b2.setSparsePoints({
        {SpeciesPair{"Fe", "Fe"}, vector{2.0}}
    });

    auto fit2 = InRamJgapFit(
        {make_shared<TwoBodyDescriptor>(desc2b2)},
        {},
        make_shared<SimpleSigmaRules>(10000000, 1, 5, 5),
        1e-8
        );

    twoAtoms.energy = 0;
    twoAtoms.atoms[0].force = Vector3{1.0, 0.0, 0.0};
    twoAtoms.atoms[1].force = Vector3{0.0, 0.0, 0.0};
    auto pot2 = fit2.fitGAP(vector{twoAtoms});

    auto coeffs2 = pot2->serialize()["descriptors"][0]["data"]["Fe,Fe"]["coefficients"];
    cout << coeffs2.dump() << endl;
    ASSERT_NEAR(coeffs2[0].get<double>() , -0.30764655994411466, 1e-6);

    // ---------------------------
    auto desc2b3 = TwoBodyDescriptor(params2b);
    desc2b3.setSparsePoints({
        {SpeciesPair{"Fe", "Fe"}, vector{2.0}}
    });

    auto fit3 = InRamJgapFit(
        {make_shared<TwoBodyDescriptor>(desc2b3)},
        {},
        make_shared<SimpleSigmaRules>(10000000, 1, 5, 5),
        1e-8
        );

    twoAtoms.energy = 0;
    twoAtoms.atoms[0].force = Vector3{1.0, 0.0, 0.0};
    twoAtoms.atoms[1].force = Vector3{0.0, 0.0, 0.0};
    twoAtoms.atoms[1].position = Vector3{1.5, 2.598, 0.0};
    auto pot3 = fit3.fitGAP(vector{twoAtoms});

    auto coeffs3 = pot3->serialize()["descriptors"][0]["data"]["Fe,Fe"]["coefficients"];
    cout << coeffs3.dump() << endl;
    ASSERT_NEAR(coeffs3[0].get<double>() , -.15382666452610533, 1e-6);

    // ---------------------------
    auto desc2b4 = TwoBodyDescriptor(params2b);
    desc2b4.setSparsePoints({
        {SpeciesPair{"Fe", "Fe"}, vector{2.0}}
    });

    auto fit4 = InRamJgapFit(
        {make_shared<TwoBodyDescriptor>(desc2b4)},
        {},
        make_shared<SimpleSigmaRules>(10000000, 2, 5, 5),
        1e-8
        );

    twoAtoms.energy = 0;
    twoAtoms.atoms[0].force = Vector3{1.0, 1.0, 1.0};
    twoAtoms.atoms[1].force = Vector3{-1.0, -1.0, -1.0};
    twoAtoms.atoms[1].position = Vector3{1.5, 2.598, 0.0};
    auto pot4 = fit4.fitGAP(vector{twoAtoms});

    auto coeffs4 = pot4->serialize()["descriptors"][0]["data"]["Fe,Fe"]["coefficients"];
    cout << coeffs4.dump() << endl;
    ASSERT_NEAR(coeffs4[0].get<double>() , -.47733536744854282, 1e-6);

    // ---------------------------
    auto desc2b5 = TwoBodyDescriptor(params2b);
    desc2b5.setSparsePoints({
        {SpeciesPair{"Fe", "Fe"}, vector{2.0, 2.5, 4.0}}
    });

    auto fit5 = InRamJgapFit(
        {make_shared<TwoBodyDescriptor>(desc2b5)},
        {},
        make_shared<SimpleSigmaRules>(1, 2, 5, 5),
        1e-8
        );

    twoAtoms.energy = 0;
    twoAtoms.atoms[0].force = Vector3{1.0, 1.0, 1.0};
    twoAtoms.atoms[1].force = Vector3{-1.0, -1.0, -1.0};
    twoAtoms.atoms[1].position = Vector3{1.5, 2.598, 0.0};
    auto pot5 = fit5.fitGAP(vector{twoAtoms});

    auto coeffs5 = pot5->serialize()["descriptors"][0]["data"]["Fe,Fe"]["coefficients"];
    cout << coeffs5.dump() << endl;
    ASSERT_NEAR(coeffs5[0].get<double>() , -.24567921601436607, 1e-6);
    ASSERT_NEAR(coeffs5[1].get<double>() , -.10998748761736651, 1e-6);
    ASSERT_NEAR(coeffs5[2].get<double>() , .38697171383300527, 1e-6);
}

TEST(TestInRamJgap, twoBodyQuipCompatibilityRealBox) {
    auto box = readXyz("test/resources/xyz-samples/FeOnly.xyz")[0];

    const auto params2b = TwoBodyDescriptorParams{
        .cutoff = 5,
        .cutoffTransitionWidth = 1,
        .kernelType = TwoBodyDescriptorParams::KernelType::GAUSS,
        .sparsificationMethod = TwoBodyDescriptorParams::SparsificationMethod::FULL_GRID_UNIFORM,
        .energyScale = 10.0,
        .lengthScale = 1.0,
    };

    auto desc2b = TwoBodyDescriptor(params2b);
    desc2b.setSparsePoints({
        {SpeciesPair{"Fe", "Fe"}, vector{2.0, 2.5, 4.0}}
    });

    auto fit = InRamJgapFit(
        {make_shared<TwoBodyDescriptor>(desc2b)},
        {},
        make_shared<SimpleSigmaRules>(0.001, 0.05, 1, 1),
        1e-8
        );
    auto pot = fit.fitGAP(vector{box});
    auto coeffs = pot->serialize()["descriptors"][0]["data"]["Fe,Fe"]["coefficients"];
    cout << coeffs.dump() << endl;

    ASSERT_NEAR(coeffs[0].get<double>() , -.91130660366840928E-002, 1e-6);
    ASSERT_NEAR(coeffs[1].get<double>() , -.30763954226713420E-003, 1e-6);
    ASSERT_NEAR(coeffs[2].get<double>() , .20225136010646777E-002, 1e-6);
}


TEST(TestInRamJgap, twoAtomsEamQuipCompatibility) {
    auto twoAtoms = AtomicStructure{
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
            }
        }
    };
    NeighbourFinder::findNeighbours(twoAtoms, 5.0);

    ////////// -----------------------------------

    constexpr auto densityParams = EamDensityCalculationParams{
        .pairFunction = EamDensityCalculationParams::PairFunction::POLYCUTOFF,
        .cutoff = 5,
        .rmin = 0,
        .speciesUse = EamDensityCalculationParams::SpeciesUse::BLIND
    };

    const auto paramsEAM = EamDescriptorParams{
        .kernelType = EamDescriptorParams::KernelType::GAUSS,
        .sparsificationMethod = EamDescriptorParams::SparsificationMethod::FULL_GRID_UNIFORM,
        .nSparsePoints = 1,
        .energyScale = 1.0,
        .lengthScale = 1.0,
        .defaultDensityCalculationParams = densityParams
    };

    shared_ptr<EamDescriptor> descEAM = make_shared<EamDescriptor>(paramsEAM);
    //descEAM->setSparsePoints({twoAtoms});
    descEAM->setSparsePoints(vector{1.0}); // \rho(3) = 0.317

    twoAtoms.energy = 1;
    twoAtoms.atoms[0].force = {1, 0, 0};
    twoAtoms.atoms[1].force = {-1, 0, 0};

    auto fit = InRamJgapFit(
        {descEAM},
        {},
        make_shared<SimpleSigmaRules>(1, 1, 1, 1),
        1e-8
        );
    auto pot = fit.fitGAP(vector{twoAtoms});
    auto coeffs = pot->serialize()["descriptors"][0]["data"]["coefficients"];
    cout << coeffs.dump() << endl;
    ASSERT_NEAR(coeffs[0].get<double>(), .17637586136984833E-001, 1e-6);
}

TEST(TestInRamJgap, eamQuipCompatibilityRealBox) {
    auto box = readXyz("test/resources/xyz-samples/FeOnly.xyz")[15];

    constexpr auto densityParams = EamDensityCalculationParams{
        .pairFunction = EamDensityCalculationParams::PairFunction::POLYCUTOFF,
        .cutoff = 5,
        .rmin = 0,
        .speciesUse = EamDensityCalculationParams::SpeciesUse::BLIND
    };

    const auto paramsEAM = EamDescriptorParams{
        .kernelType = EamDescriptorParams::KernelType::GAUSS,
        .sparsificationMethod = EamDescriptorParams::SparsificationMethod::FULL_GRID_UNIFORM,
        .nSparsePoints = 3,
        .energyScale = 2.0,
        .lengthScale = 3.0,
        .defaultDensityCalculationParams = densityParams
    };

    shared_ptr<EamDescriptor> descEAM = make_shared<EamDescriptor>(paramsEAM);
    descEAM->setSparsePoints(vector{1.0, 2.5, 4.0});

    auto fit = InRamJgapFit(
        {descEAM},
        {},
        make_shared<SimpleSigmaRules>(0.001, 0.05, 1, 1),
        1e-8
        );

    auto pot = fit.fitGAP(vector{box});
    auto coeffs = pot->serialize()["descriptors"][0]["data"]["coefficients"];
    cout << coeffs.dump() << endl;
/*
    <sparseX i="1" alpha="-5.0619173053329485" sparseCutoff="1.0000000000000000"/>
    <sparseX i="2" alpha="-4.4169841885357615" sparseCutoff="1.0000000000000000"/>
    <sparseX i="3" alpha="-.16117464517539212" sparseCutoff="1.0000000000000000"/>
    </gpCoordinates>
*/
    ASSERT_NEAR(coeffs[0].get<double>(), -5.0619173053329485, 1e-6);
    ASSERT_NEAR(coeffs[1].get<double>(), -4.4169841885357615, 1e-6);
    ASSERT_NEAR(coeffs[2].get<double>(), -.16117464517539212, 1e-6);
}

TEST(TestInRamJgap, equilateralTriangle3bQuipCompatibility) {
    setupEquilateralTriangle();

    equilateralTriangle.energy = 1.0;
    equilateralTriangle.atoms[0].force = {1, 0, 0};
    equilateralTriangle.atoms[1].force = {0, -1, 0};
    equilateralTriangle.atoms[2].force = {0.5, 0.5, 0};

    const auto params3b = ThreeBodyDescriptorParams{
        .cutoff = 10.0,
        .cutoffTransitionWidth = 0.6,
        .kernelType = ThreeBodyDescriptorParams::KernelType::GAUSS,
        .sparsificationMethod = ThreeBodyDescriptorParams::SparsificationMethod::FULL_GRID_UNIFORM,
        .energyScale = 1.0,
        .lengthScale = 1.0
    };

    auto desc3b = make_shared<ThreeBodyDescriptor>(params3b);
    desc3b->setSparsePoints(map<SpeciesTriplet, vector<Vector3>>{
        {
            SpeciesTriplet{.root = "Fe", .nodes = SpeciesPair{"Fe", "Fe"}},
            {Vector3{6, 0, 3}}
        }
    });

    auto fit = InRamJgapFit(
    {desc3b},
    {},
    make_shared<SimpleSigmaRules>(1, 1, 1, 1),
    1e-8
    );

    auto pot = fit.fitGAP(vector{equilateralTriangle});
    // cout << pot->serialize();
    auto coeffs = pot->serialize()["descriptors"][0]["data"]["Fe,Fe,Fe"]["coefficients"];
    cout << coeffs.dump() << endl;
    ASSERT_NEAR(coeffs[0].get<double>(), .15381640115291151, 1e-8);
}

AtomicStructure pythagorian;

void initPythagorian() {
    pythagorian = AtomicStructure{
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
    NeighbourFinder::findNeighbours(pythagorian, 10.0);
}

TEST(TestInRamJgap, pythagorian3bQuipCompatibility) {
    initPythagorian();

    pythagorian.energy = 1.0;
    pythagorian.atoms[0].force = {1, 0, 0};
    pythagorian.atoms[1].force = {0, -1, 0};
    pythagorian.atoms[2].force = {0.5, 0.5, 0};

    const auto params3b = ThreeBodyDescriptorParams{
        .cutoff = 10.0,
        .cutoffTransitionWidth = 0.6,
        .kernelType = ThreeBodyDescriptorParams::KernelType::GAUSS,
        .sparsificationMethod = ThreeBodyDescriptorParams::SparsificationMethod::FULL_GRID_UNIFORM,
        .energyScale = 1.0,
        .lengthScale = 1.0
    };

    auto desc3b = make_shared<ThreeBodyDescriptor>(params3b);
    desc3b->setSparsePoints(map<SpeciesTriplet, vector<Vector3>>{
        {
            SpeciesTriplet{.root = "Fe", .nodes = SpeciesPair{"Fe", "Fe"}},
            {Vector3{6, 0, 3}}
        }
    });

    auto fit = InRamJgapFit(
    {desc3b},
    {},
    make_shared<SimpleSigmaRules>(1, 1, 1, 1),
    1e-8
    );

    auto pot = fit.fitGAP(vector{pythagorian});
    // cout << pot->serialize();
    auto coeffs = pot->serialize()["descriptors"][0]["data"]["Fe,Fe,Fe"]["coefficients"];
    cout << coeffs.dump() << endl;
    ASSERT_NEAR(coeffs[0].get<double>(), -.13044210897182607, 1e-8);
}

TEST(TestInRamJgap, hard3bQuipCompatibility) {
    auto box = readXyz("test/resources/xyz-samples/FeOnly.xyz")[15];

    const auto params3b = ThreeBodyDescriptorParams{
        .cutoff = 4.5,
        .cutoffTransitionWidth = 0.6,
        .kernelType = ThreeBodyDescriptorParams::KernelType::GAUSS,
        .sparsificationMethod = ThreeBodyDescriptorParams::SparsificationMethod::FULL_GRID_UNIFORM,
        .energyScale = 3.0,
        .lengthScale = 0.5
    };

    auto desc3b = make_shared<ThreeBodyDescriptor>(params3b);
    desc3b->setSparsePoints(map<SpeciesTriplet, vector<Vector3>>{
        {
            SpeciesTriplet{.root = "Fe", .nodes = SpeciesPair{"Fe", "Fe"}},
            {
                Vector3{6, 0, 3},
                Vector3{8, 2, 1},
                Vector3{5, 3, 8}
            }
        }
    });

    auto fit = InRamJgapFit(
    {desc3b},
    {},
    make_shared<SimpleSigmaRules>(0.001, 0.05, 1, 1),
    1e-8
    );

    auto pot = fit.fitGAP(vector{box});
    // cout << pot->serialize();
    auto coeffs = pot->serialize()["descriptors"][0]["data"]["Fe,Fe,Fe"]["coefficients"];
    cout << coeffs.dump() << endl;
    /*      <sparseX i="1" alpha="-.63024066682173821E-001" sparseCutoff="1.0000000000000000"/>
      <sparseX i="2" alpha="20.858471601167693" sparseCutoff="1.0000000000000000"/>
      <sparseX i="3" alpha="-190.02399797873642" sparseCutoff="1.0000000000000000"/>*/
    ASSERT_NEAR(coeffs[0].get<double>(), -.63024066682173821E-001, 1e-8);
    ASSERT_NEAR(coeffs[1].get<double>(), 20.858471601167693, 1e-8);
    ASSERT_NEAR(coeffs[2].get<double>(), -190.02399797873642, 1e-8);
}