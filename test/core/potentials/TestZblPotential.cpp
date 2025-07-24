#include <gtest/gtest.h>
#include <core/potentials/ZblPotential.hpp>

#include <fstream>
#include <iostream>

#include "core/neighbours/NeighbourFinder.hpp"

using namespace std;
using namespace jgap;

TEST(TestZblPotential, DimersData) {
    auto zblPot = ZblPotential(nlohmann::json::parse("{}"));

    string xyzDataFn = "test/resources/zbl/quip-test.out.xyz";
    auto structs = readXyz(xyzDataFn);
    NeighbourFinder::findNeighbours(structs, 5);

    vector selected = {structs[1]};
    for (int i: {24, 25, 26, 28}) {
        for (int j: {24, 25, 26, 28}) {
            if (j < i) {
                continue;
            }
            ofstream ffout(format("jgap-dimers-{}-{}.dat",i, j));
            for (double d = 1; d <= 2; d+=0.01) {
                structs[1].positions[1].x = d;
                structs[1].species[0] = Z_inverse[i];
                structs[1].species[1] = Z_inverse[j];
                NeighbourFinder::findNeighbours(selected, 3);

                auto pred = zblPot.predict(structs[1]).energy.value();
                ffout <<d << " " << pred << "\n";
            }
            ffout.close();
        }
    }
}

// TODO: re-fit/check quip cutoff
void testVsQuipPairpot(const string& xyzDataFn) {

    auto zblPot = ZblPotential(nlohmann::json::parse("{}"));

    auto structs = readXyz(xyzDataFn);
    NeighbourFinder::findNeighbours(structs, 5);

    for (size_t i = 0; i < structs.size(); i++) {

        auto pred = zblPot.predict(structs[i]);

        if (abs(structs[i].energy.value()) < 1.0) {
            cout << i << endl;
            ASSERT_NEAR(pred.energy.value(), structs[i].energy.value(), 0.08);
        } else {
            ASSERT_NEAR((pred.energy.value()-structs[i].energy.value())/structs[i].energy.value(), 0, 0.08);
        }

        for (size_t j = 0; j < structs[i].size(); j++) {
            double closestNeighbour = 5;
            for (const NeighbourData& neigh: structs[i][j].neighbours()) {
                if (neigh.distance < closestNeighbour) {
                    closestNeighbour = neigh.distance;
                }
            }

            Vector3 fDiff = structs[i][j].force() - pred.forces.value()[j];
            if (closestNeighbour < 1.5) {
                ASSERT_NEAR(fDiff.len() / structs[i][j].force().len(), 0, 0.4);
            }
        }
    }
}

TEST(TestZblPotential, CompareIter33Test) {
    auto filename = "test/resources/zbl/quip-test.out.xyz";
    testVsQuipPairpot(filename);
}

TEST(TestZblPotential, CompareIter33Train) {
    auto filename = "test/resources/zbl/quip-train.out.xyz";
    testVsQuipPairpot(filename);
}