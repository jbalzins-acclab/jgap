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
#include <execinfo.h>

void print_backtrace() {
    void* callstack[128];
    int frames = backtrace(callstack, 128);
    char** symbols = backtrace_symbols(callstack, frames);
    std::cerr << "Stack trace:\n";
    for (int i = 0; i < frames; ++i) {
        std::cerr << symbols[i] << "\n";
    }
    free(symbols);
}

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

    } catch (exception& e) {
        jgap::CurrentLogger::get()->error("Fail: " + string(e.what()));
        print_backtrace();
        throw;
    }

    return EXIT_SUCCESS;
}
