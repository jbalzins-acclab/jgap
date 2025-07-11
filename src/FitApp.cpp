#include <core/fit/InRamJgapFit.hpp>
#include <utils/Utils.hpp>
#include <iostream>
#include <fstream>
#include <nlohmann/json.hpp>
#include <format>
#include <cmath>
#include <thread>
#include <Eigen/Dense>
#include <ParserRegistryAuto.hpp>
#include <tbb/parallel_for_each.h>

using namespace std;

int main(int argc, char** argv) {

    Eigen::setNbThreads(thread::hardware_concurrency()); // NOLINT(*-narrowing-conversions)

    if (argc != 2) {
        jgap::CurrentLogger::get()->error("Expected 1 argument: \n fit param file");
        return EXIT_FAILURE;
    }

    try {

        // ------------------------ READ PARAMS AND PREPARE -------------------------------

        string paramFileName = argv[1];
        ifstream paramFile(paramFileName);
        if (!paramFile.is_open()) {
            jgap::CurrentLogger::get()->error(format("Cannot open param file {}", paramFileName));
            return EXIT_FAILURE;
        }
        nlohmann::json fitParams;
        paramFile >> fitParams;

        string outputFileName = fitParams["output_file"];
        ofstream outFileTest(outputFileName);
        if (!outFileTest.is_open()) {
            jgap::CurrentLogger::get()->error(format("Cannot open output file {}", outputFileName), true);
            return EXIT_FAILURE;
        }
        outFileTest.close();

        fitParams["type"] = fitParams.value("type", "composite");
        auto fit = jgap::ParserRegistry<jgap::Fit>::get(fitParams);

        auto trainingData = jgap::readXyz(fitParams["training_data_xyz"]);

        // ------------------------ FIT -------------------------------

        jgap::CurrentLogger::get()->info(format(
            "Fitting {} potential with params from file {}: {}",
            fitParams["type"].dump(), paramFileName, fitParams.dump()
            ));
        shared_ptr<jgap::Potential> resultingPotential = fit->fit(trainingData);
        jgap::CurrentLogger::get()->info("Main fit done");

        // ------------------------ SAVE -------------------------------

        jgap::CurrentLogger::get()->info(format("Saving resulting potential data to {}", outputFileName));
        ofstream outFile(outputFileName);
        outFile << resultingPotential->serialize().dump(4) << endl;
        outFile.flush();
        outFile.close();

        // ------------------------ TEST IF REQUESTED -------------------------------

        if (!fitParams.contains("testing_data_xyzs")) return EXIT_SUCCESS;

        for (const auto& testFileName: fitParams["testing_data_xyzs"]) {

            jgap::CurrentLogger::get()->info(format("Testing at {}", testFileName.get<string>()));

            auto testingData = jgap::readXyz(testFileName);
            jgap::NeighbourFinder::findNeighbours(testingData, resultingPotential->getCutoff());

            vector<double> dE, dF;

            for (const auto& structure : testingData) {
                if (dE.size() % max(testingData.size()/100, 1uz) == 0) {
                    jgap::CurrentLogger::get()->debug(format(
                        "Testing at {}%",
                        dE.size() * 100 / testingData.size()
                        ));
                }

                auto prediction = resultingPotential->predict(structure);
                if (structure.energy.has_value() && prediction.energy.has_value()) {
                    dE.push_back((structure.energy.value() - prediction.energy.value())
                                        / static_cast<double>(structure.atoms.size()));
                }

                if (!prediction.forces.has_value()) continue;

                for (size_t i = 0; i < structure.atoms.size(); i++) {
                    if (structure.atoms[i].force.has_value()) {
                        dF.push_back((structure.atoms[i].force.value() - prediction.forces.value()[i]).norm());
                    }
                }
            }

            jgap::CurrentLogger::get()->info(format("E_rmse = {}eV/A", jgap::rms(dE)));
            jgap::CurrentLogger::get()->info(format("E_std = {}eV", jgap::std(dE)));
            jgap::CurrentLogger::get()->info(format("F_rmse = {}eV/A", jgap::rms(dF)));
            jgap::CurrentLogger::get()->info(format("F_std = {}eV/A", jgap::std(dF)));
        }

    } catch (exception& e) {
        jgap::CurrentLogger::get()->error("Fail: " + string(e.what()));
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
