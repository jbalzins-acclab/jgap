#include <core/fit/QRGapFit.hpp>
#include <utils/Utils.hpp>
#include <iostream>
#include <fstream>
#include <nlohmann/json.hpp>
#include <format>
#include <ParserRegistryAuto.hpp>
#include <Version.hpp>
#include <tbb/parallel_for.h>

using namespace std;

int main(int argc, char** argv) {

    jgap::CurrentLogger::get()->info("jGAP predict v{}", JGAP_VERSION);

    if (argc != 4) {
        jgap::CurrentLogger::get()->error(
            "Expected 3 arguments: \n 1. pot-file(.jgap.json) \n 2. to-be-predicted.xyz \n 3. output.xyz"
            );
        return EXIT_FAILURE;
    }

    try {

        // ------------------------ READ PARAMS AND PREPARE -------------------------------
        jgap::CurrentLogger::get()->info("Using {} potential to predict on data from {}", argv[1], argv[2]);

        string potentialFileName = argv[1];
        ifstream potentialParamFile(potentialFileName);
        if (!potentialParamFile.is_open()) {
            jgap::CurrentLogger::get()->error("Cannot open pot-param file {}", potentialFileName);
            return EXIT_FAILURE;
        }
        nlohmann::json potParams;
        potentialParamFile >> potParams;
        if (!potParams.contains("type")) {
            potParams["type"] = "composite";
        }

        jgap::CurrentLogger::get()->info("Loading potential");
        auto potential = jgap::ParserRegistry<jgap::Potential>::get(potParams);

        jgap::CurrentLogger::get()->info("Loading data");
        auto toBePredicted = jgap::readXyz(argv[2]);

        // ------------------------ PREDICT IN PARALLEL -------------------------------
        jgap::CurrentLogger::get()->info("Evaluation started");
        jgap::NeighbourFinder::findNeighbours(toBePredicted, potential->getCutoff());
        vector<jgap::PotentialPrediction> predictions(toBePredicted.size());
        tbb::parallel_for(0uz, toBePredicted.size(), [&](const size_t i) {
            predictions[i] = potential->predict(toBePredicted[i]);
        });
        jgap::CurrentLogger::get()->info("Evaluation finished");

        // ------------------------ GATHER PREDICTIONS AND SAVE -------------------------------
        vector<jgap::AtomicStructure> result;
        for (size_t i = 0; i < toBePredicted.size(); i++) {
            auto structure = toBePredicted[i];
            structure.setEnergyData(predictions[i]);
            result.push_back(structure);
        }

        jgap::writeXyz(argv[3], result);
        jgap::CurrentLogger::get()->info("Finished predicting");

    } catch (exception& e) {
        jgap::CurrentLogger::get()->error("Fail: {}", e.what());
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
