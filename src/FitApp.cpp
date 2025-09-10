#include <core/fit/InRamJgapFit.hpp>
#include <utils/Utils.hpp>
#include <fstream>
#include <nlohmann/json.hpp>
#include <thread>
#include <Eigen/Dense>
#include <ParserRegistryAuto.hpp>
#include <tbb/parallel_for_each.h>

#include <Version.hpp>

using namespace std;

int main(int argc, char** argv) {
    jgap::CurrentLogger::get()->info("jGAP fit v{}", JGAP_VERSION);

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
            jgap::CurrentLogger::get()->error("Cannot open param file {}", paramFileName);
            return EXIT_FAILURE;
        }
        nlohmann::json fitParams;
        paramFile >> fitParams;

        string outputFileName = fitParams["output_file"];
        jgap::CurrentLogger::get()->info("Output file name: {}", outputFileName);

        jgap::CurrentLogger::get()->info("Picking fit logic");
        fitParams["type"] = fitParams.value("type", "composite");
        auto fit = jgap::ParserRegistry<jgap::Fit>::get(fitParams);

        jgap::CurrentLogger::get()->info("Reading training data");
        auto trainingData = jgap::readXyz(fitParams["training_data_xyz"]);

        if (!fitParams.value("use_virials", true)) {
            jgap::CurrentLogger::get()->warn("Not using virials in training data");

            for (auto& structure: trainingData) {
                structure.virials.reset();
            }
        }

        // ------------------------ FIT -------------------------------

        jgap::CurrentLogger::get()->info(
            "Fitting {} potential with params from file {}: {}",
            fitParams["type"].dump(), paramFileName, fitParams.dump()
            );
        shared_ptr<jgap::Potential> resultingPotential = fit->fit(trainingData);
        jgap::CurrentLogger::get()->info("Main fit done");

        // ------------------------ SAVE -------------------------------

        jgap::CurrentLogger::get()->info("Saving resulting potential data to {}", outputFileName);
        auto output = resultingPotential->serialize();
        output["type"] = resultingPotential->getType();

        ofstream outFile(outputFileName);
        outFile << output.dump(4) << endl;
        outFile.flush();
        outFile.close();

        jgap::CurrentLogger::get()->info("Fit complete");

    } catch (exception& e) {
        jgap::CurrentLogger::get()->error("Fail: " + string(e.what()));
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
