#include <gtest/gtest.h>

#include "core/fit/InRamJgapFit.hpp"
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

    const auto params = nlohmann::json::parse(R"(
    {
        "descriptors": {
            "2b_test": {
                "type": "2b",
                "kernel": {
                    "type": "squared_exp",
                    "length_scale": 1.0,
                    "energy_scale": 1.0,
                    "cutoff": {
                        "type": "coscutoff",
                        "r_min": 9.3,
                        "cutoff": 10.0
                    }
                },
                "sparse_data": {
                    "Fe,Fe": {
                        "sparse_points": [3.0]
                    }
                }
            }
        },
        "jitter": 1e-8,
        "sigma_rules": {
            "type": "simple",
            "E_per_root_n_atoms": 3.0,
            "F_component": 10.0,
            "liquid": 5,
            "short_range": 5
        }
    }
    )");

    auto fit = InRamJgapFit(params);

    equilateralTriangle.energy = 1;
    equilateralTriangle.atoms[0].force = Vector3{0.0, 0.0, 0.0};
    equilateralTriangle.atoms[1].force = Vector3{0.0, 0.0, 0.0};
    equilateralTriangle.atoms[2].force = Vector3{0.0, 0.0, 0.0};
    //NeighbourFinder::findNeighbours(equilateralTriangle, 10.0);
    //auto a = desc2b.covariate(equilateralTriangle);
    auto pot = fit.fit(vector{equilateralTriangle});

    auto res = pot->serialize();
    auto coeffs = res["descriptors"]["2b_test"]["sparse_data"]["Fe,Fe"]["coefficients"];
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

nlohmann::json makeParams2bDesc(double theta, double delta, double rMin, double cutoff, vector<double> sparsePts) {

    return nlohmann::json{
        {"type", "2b"},
        {"kernel", {
            {"type", "squared_exp"},
            {"length_scale", theta},
            {"energy_scale", delta},
            {"cutoff", {
                    {"type", "coscutoff"},
                    {"r_min", rMin},
                    {"cutoff", cutoff}
                    }
                }
            }
        },
        {"sparse_data", {{"Fe,Fe", {{"sparse_points", sparsePts}}}}}
    };
}

nlohmann::json makeSimpleSigmaRules(double E, double F, double liquid, double shortRange) {
    return nlohmann::json{
        {"type", "simple"},
        {"E_per_root_n_atoms", E},
        {"F_component", F},
        {"liquid", liquid},
        {"short_range", shortRange}
    };
}

nlohmann::json makeFitParams(nlohmann::json desc, string descName, nlohmann::json srules) {
    return nlohmann::json{
        {"descriptors", {
                {descName, desc}
            }
        },
        {"jitter", 1e-8},
        {"sigma_rules", srules}
    };
}

AtomicStructure twoAtoms;
void setupTwoAtoms() {
    twoAtoms = AtomicStructure{
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
}

TEST(TestInRamJgap, twoAtomsWithForceQuipCompatibility1) {
    setupTwoAtoms();
    const auto params = makeFitParams(
            makeParams2bDesc(1.0, 1.0, 9.3, 10.0, vector{2.0, 4.0}),
            "2b_test",
            makeSimpleSigmaRules(1, 1, 5, 5)
        );

    auto fit = InRamJgapFit(params);

    twoAtoms.energy = 1;
    twoAtoms.atoms[0].force = Vector3{0.0, 0.0, 0.0};
    twoAtoms.atoms[1].force = Vector3{0.0, 0.0, 0.0};
    auto pot = fit.fit(vector{twoAtoms});

    auto coeffs = pot->serialize()["descriptors"]["2b_test"]["sparse_data"]["Fe,Fe"]["coefficients"];
    cout << coeffs.dump() << endl;
    auto a0 = coeffs[0].get<double>(), a1 = coeffs[1].get<double>();
    ASSERT_NEAR(a0, 0.23266775926220731, 1e-6);
    ASSERT_NEAR(a1, 0.23266775926220731, 1e-6);
}

TEST(TestInRamJgap, twoAtomsWithForceQuipCompatibility2) {

    setupTwoAtoms();
    const auto params = makeFitParams(
            makeParams2bDesc(1.0, 1.0, 9.3, 10.0, vector{2.0}),
            "2b_test",
            makeSimpleSigmaRules(10000000, 1, 5, 5)
        );

    auto fit = InRamJgapFit(params);

    twoAtoms.energy = 0;
    twoAtoms.atoms[0].force = Vector3{1.0, 0.0, 0.0};
    twoAtoms.atoms[1].force = Vector3{0.0, 0.0, 0.0};
    auto pot = fit.fit(vector{twoAtoms});

    auto coeffs = pot->serialize()["descriptors"]["2b_test"]["sparse_data"]["Fe,Fe"]["coefficients"];
    cout << coeffs.dump() << endl;
    ASSERT_NEAR(coeffs[0].get<double>() , -0.30764655994411466, 1e-6);
}

TEST(TestInRamJgap, twoAtomsWithForceQuipCompatibility3) {
    setupTwoAtoms();
    const auto params = makeFitParams(
            makeParams2bDesc(1.0, 1.0, 9.3, 10.0, vector{2.0}),
            "2b_test",
            makeSimpleSigmaRules(10000000, 1, 5, 5)
        );

    auto fit = InRamJgapFit(params);

    twoAtoms.energy = 0;
    twoAtoms.atoms[0].force = Vector3{1.0, 0.0, 0.0};
    twoAtoms.atoms[1].force = Vector3{0.0, 0.0, 0.0};
    twoAtoms.atoms[1].position = Vector3{1.5, 2.598, 0.0};
    auto pot = fit.fit(vector{twoAtoms});

    auto coeffs = pot->serialize()["descriptors"]["2b_test"]["sparse_data"]["Fe,Fe"]["coefficients"];
    cout << coeffs.dump() << endl;
    ASSERT_NEAR(coeffs[0].get<double>() , -.15382666452610533, 1e-6);
}

TEST(TestInRamJgap, twoAtomsWithForceQuipCompatibility4) {
    setupTwoAtoms();
    const auto params = makeFitParams(
            makeParams2bDesc(1.0, 1.0, 9.3, 10.0, vector{2.0}),
            "2b_test",
            makeSimpleSigmaRules(10000000, 2, 5, 5)
        );

    auto fit = InRamJgapFit(params);

    twoAtoms.energy = 0;
    twoAtoms.atoms[0].force = Vector3{1.0, 1.0, 1.0};
    twoAtoms.atoms[1].force = Vector3{-1.0, -1.0, -1.0};
    twoAtoms.atoms[1].position = Vector3{1.5, 2.598, 0.0};
    auto pot = fit.fit(vector{twoAtoms});

    auto coeffs = pot->serialize()["descriptors"]["2b_test"]["sparse_data"]["Fe,Fe"]["coefficients"];
    cout << coeffs.dump() << endl;
    ASSERT_NEAR(coeffs[0].get<double>() , -.47733536744854282, 1e-6);
}

TEST(TestInRamJgap, twoAtomsWithForceQuipCompatibility5) {
    setupTwoAtoms();
    const auto params = makeFitParams(
            makeParams2bDesc(1.0, 1.0, 9.3, 10.0, vector{2.0, 2.5, 4.0}),
            "2b_test",
            makeSimpleSigmaRules(1, 2, 5, 5)
        );

    auto fit = InRamJgapFit(params);

    twoAtoms.energy = 0;
    twoAtoms.atoms[0].force = Vector3{1.0, 1.0, 1.0};
    twoAtoms.atoms[1].force = Vector3{-1.0, -1.0, -1.0};
    twoAtoms.atoms[1].position = Vector3{1.5, 2.598, 0.0};
    auto pot = fit.fit(vector{twoAtoms});

    auto coeffs = pot->serialize()["descriptors"]["2b_test"]["sparse_data"]["Fe,Fe"]["coefficients"];
    cout << coeffs.dump() << endl;
    ASSERT_NEAR(coeffs[0].get<double>() , -.24567921601436607, 1e-6);
    ASSERT_NEAR(coeffs[1].get<double>() , -.10998748761736651, 1e-6);
    ASSERT_NEAR(coeffs[2].get<double>() , .38697171383300527, 1e-6);
}

TEST(TestInRamJgap, twoBodyQuipCompatibilityRealBox) {
    auto box = readXyz("test/resources/xyz-samples/FeOnly.xyz")[0];
    const auto params = makeFitParams(
            makeParams2bDesc(1.0, 10.0, 4.0, 5.0, vector{2.0, 2.5, 4.0}),
            "2b_test",
            makeSimpleSigmaRules(0.001, 0.05, 1, 1)
        );

    auto fit = InRamJgapFit(params);

    auto pot = fit.fit(vector{box});
    auto coeffs = pot->serialize()["descriptors"]["2b_test"]["sparse_data"]["Fe,Fe"]["coefficients"];
    cout << coeffs.dump() << endl;

    ASSERT_NEAR(coeffs[0].get<double>() , -.91130660366840928E-002, 1e-6);
    ASSERT_NEAR(coeffs[1].get<double>() , -.30763954226713420E-003, 1e-6);
    ASSERT_NEAR(coeffs[2].get<double>() , .20225136010646777E-002, 1e-6);
}

nlohmann::json makeParamsEamDesc(double theta, double delta, double rMin, double cutoff, vector<double> sparsePts) {

    return nlohmann::json{
        {"type", "eam"},
        {"kernel", {
            {"type", "squared_exp"},
            {"length_scale", theta},
            {"energy_scale", delta}
            }
        },
        {"sparse_data", {{"Fe", {{"sparse_points", sparsePts}}}}},
        {"pair_functions", {
                {
                    {"type", "polycutoff"},
                    {"cutoff", cutoff},
                    {"r_min", rMin}
                }
            }
        }
    };
}

TEST(TestInRamJgap, twoAtomsEamQuipCompatibility) {
    setupTwoAtoms();

    twoAtoms.energy = 1;
    twoAtoms.atoms[0].force = {1, 0, 0};
    twoAtoms.atoms[1].force = {-1, 0, 0};

    const auto params = makeFitParams(
            makeParamsEamDesc(1.0, 1.0, 0.0, 5.0, vector{1.0}),
            "eam_test",
            makeSimpleSigmaRules(1, 1, 1, 1)
        );
    auto fit = InRamJgapFit(params);

    auto pot = fit.fit(vector{twoAtoms});
    auto coeffs = pot->serialize()["descriptors"]["eam_test"]["sparse_data"]["Fe"]["coefficients"];
    cout << coeffs.dump() << endl;
    ASSERT_NEAR(coeffs[0].get<double>(), .17637586136984833E-001, 1e-6);
}

TEST(TestInRamJgap, eamQuipCompatibilityRealBox) {
    auto box = readXyz("test/resources/xyz-samples/FeOnly.xyz")[15];

    const auto params = makeFitParams(
        makeParamsEamDesc(3.0, 2.0, 0.0, 5.0, vector{1.0, 2.5, 4.0}),
        "eam_test",
        makeSimpleSigmaRules(0.001, 0.05, 1, 1)
    );
    auto fit = InRamJgapFit(params);

    auto pot = fit.fit(vector{box});
    auto coeffs = pot->serialize()["descriptors"]["eam_test"]["sparse_data"]["Fe"]["coefficients"];
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

nlohmann::json makeParams3bDesc(double theta, double delta, double rMin, double cutoff, vector<Vector3> sparsePts) {

    auto spConv = nlohmann::json::array();
    for (const auto &sp: sparsePts) {
        spConv.push_back({
            {"x", sp.x},
            {"y", sp.y},
            {"z", sp.z}
        });
    }

    return nlohmann::json{
            {"type", "3b"},
            {"kernel", {
                {"type", "squared_exp"},
                {"length_scale", theta},
                {"energy_scale", delta},
                {"cutoff", {
                        {"type", "coscutoff"},
                        {"r_min", rMin},
                        {"cutoff", cutoff}
                        }
                    }
                }
            },
            {"sparse_data", {{"Fe,Fe,Fe", {{"sparse_points", spConv}}}}}
    };
}

TEST(TestInRamJgap, equilateralTriangle3bQuipCompatibility) {
    setupEquilateralTriangle();

    equilateralTriangle.energy = 1.0;
    equilateralTriangle.atoms[0].force = {1, 0, 0};
    equilateralTriangle.atoms[1].force = {0, -1, 0};
    equilateralTriangle.atoms[2].force = {0.5, 0.5, 0};

    const auto params = makeFitParams(
            makeParams3bDesc(1.0, 1.0, 9.4, 10.0, vector{Vector3{6.0, 0.0, 3.0}}),
            "3b_test",
            makeSimpleSigmaRules(1, 1, 1, 1)
        );

    auto fit = InRamJgapFit(params);

    auto pot = fit.fit(vector{equilateralTriangle});
    // cout << pot->serialize();
    auto coeffs = pot->serialize()["descriptors"]["3b_test"]["sparse_data"]["Fe,Fe,Fe"]["coefficients"];
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

    const auto params = makeFitParams(
            makeParams3bDesc(1.0, 1.0, 9.4, 10.0, vector{Vector3{6.0, 0.0, 3.0}}),
            "3b_test",
            makeSimpleSigmaRules(1, 1, 1, 1)
        );

    auto fit = InRamJgapFit(params);

    auto pot = fit.fit(vector{pythagorian});
    // cout << pot->serialize();
    auto coeffs = pot->serialize()["descriptors"]["3b_test"]["sparse_data"]["Fe,Fe,Fe"]["coefficients"];
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

    const auto params = makeFitParams(
    makeParams3bDesc(0.5, 3.0, 3.9, 4.5, vector{
                Vector3{6, 0, 3},
                Vector3{8, 2, 1},
                Vector3{5, 3, 8}
            }),
            "3b_test",
            makeSimpleSigmaRules(0.001, 0.05, 1, 1)
        );

    auto fit = InRamJgapFit(params);

    auto pot = fit.fit(vector{box});
    // cout << pot->serialize();
    auto coeffs = pot->serialize()["descriptors"]["3b_test"]["sparse_data"]["Fe,Fe,Fe"]["coefficients"];
    cout << coeffs.dump() << endl;
    /*      <sparseX i="1" alpha="-.63024066682173821E-001" sparseCutoff="1.0000000000000000"/>
      <sparseX i="2" alpha="20.858471601167693" sparseCutoff="1.0000000000000000"/>
      <sparseX i="3" alpha="-190.02399797873642" sparseCutoff="1.0000000000000000"/>*/
    ASSERT_NEAR(coeffs[0].get<double>(), -.63024066682173821E-001, 1e-8);
    ASSERT_NEAR(coeffs[1].get<double>(), 20.858471601167693, 1e-8);
    ASSERT_NEAR(coeffs[2].get<double>(), -190.02399797873642, 1e-8);
}