#include <gtest/gtest.h>
#include <core/potentials/ZblPotential.hpp>

#include <fstream>
#include <iostream>

#include "core/neighbours/NeighbourFinder.hpp"

using namespace std;
using namespace jgap;

TEST(TestZblPotential, DimersData) {
    ifstream file("resources/dmol-screening-fit/dmol-fit.json");
    nlohmann::json dmolData;
    file >> dmolData;
    auto zblPot = ZblPotential(dmolData);

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
                structs[1].atoms[1].position.x = d;
                structs[1].atoms[0].species = Z_inverse[i];
                structs[1].atoms[1].species = Z_inverse[j];
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

    ifstream file("resources/dmol-screening-fit/dmol-fit.json");
    nlohmann::json dmolData;
    file >> dmolData;
    auto zblPot = ZblPotential(dmolData);

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

        for (size_t j = 0; j < structs[i].atoms.size(); j++) {
            double closestNeighbour = 5;
            for (NeighbourData& neigh: structs[i].atoms[j].neighbours.value()) {
                if (neigh.distance < closestNeighbour) {
                    closestNeighbour = neigh.distance;
                }
            }

            Vector3 fDiff = structs[i].atoms[j].force.value() - pred.forces.value()[j];
            if (closestNeighbour < 1.5) {
                ASSERT_NEAR(fDiff.norm() / structs[i].atoms[j].force.value().norm(), 0, 0.4);
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