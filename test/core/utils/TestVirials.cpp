#include <gtest/gtest.h>
#include <string>
#include <fstream>
#include <nlohmann/json.hpp>
#include <ParserRegistryAuto.hpp>

#include "core/descriptors/TwoBodyDescriptor.hpp"
#include "utils/Utils.hpp"

using namespace std;
using namespace jgap;

TEST(TestVirials, quip) {

    string xyzDataFn = "test/resources/xyz-samples/fe-only-compact.xyz";
    auto box = readXyz(xyzDataFn, 10.0)[0];

    ifstream fIn("resources/jgap-pots/3b-desc.json");
    nlohmann::json potParams;
    fIn >> potParams;

    CompositePotential desc(potParams);

    auto res = desc.predict(box);
    0;
}
