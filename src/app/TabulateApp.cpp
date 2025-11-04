#include <core/fit/QRGapFit.hpp>
#include <utils/Utils.hpp>
#include <iostream>
#include <fstream>
#include <nlohmann/json.hpp>
#include <format>
#include <ParserRegistryAuto.hpp>
#include <Version.hpp>
#include <tbb/parallel_for.h>

#include "../../include/core/tabulation/TabGapFit.hpp"
#include "io/convert/QuipXmlConverter.hpp"

using namespace std;

int main(int argc, char** argv) {

    jgap::CurrentLogger::get()->info("jGAP tabulate v{}", JGAP_VERSION);

    if (argc != 2) {
        jgap::CurrentLogger::get()->error(
            "Expected 1 arguments: \n 1. <<tabulation-params>>.json"
            );
        return EXIT_FAILURE;
    }

    try {

        // ------------------------ READ PARAMS AND PREPARE -------------------------------
        jgap::CurrentLogger::get()->info("Tabulation as specified in: {}", argv[1]);

        string paramFileName = argv[1];
        ifstream paramFile(paramFileName);
        if (!paramFile.is_open()) {
            jgap::CurrentLogger::get()->error("Cannot open tabulation-param file {}", paramFileName);
            return EXIT_FAILURE;
        }

        nlohmann::json tabulationParams;
        paramFile >> tabulationParams;

        const string outputFilePrefix = tabulationParams["output_file_prefix"];
        const string potentialFileName = tabulationParams["potential_file"];

        shared_ptr<jgap::Potential> potential;
        if (potentialFileName.ends_with(".xml")) {

            pugi::xml_document quipDocument;
            if (!quipDocument.load_file(potentialFileName.c_str())) {
                jgap::CurrentLogger::get()->error("Cannot input quip file: {}",  potentialFileName);
                return EXIT_FAILURE;
            }

            jgap::CurrentLogger::get()->info("Converting from QUIP xml");
            potential = jgap::QuipXmlConverter::transform(quipDocument.document_element());

        } else {
            ifstream potentialParamFile(potentialFileName);
            if (!potentialParamFile.is_open()) {
                jgap::CurrentLogger::get()->error("Cannot open pot-param file {}", potentialFileName);
                return EXIT_FAILURE;
            }
            nlohmann::json potParams;
            potentialParamFile >> potParams;
            potential = jgap::ParserRegistry<jgap::Potential>::get(potParams);
        }

        jgap::CurrentLogger::get()->info("Tabulating potential: {}", potential->serialize().dump());
        jgap::CurrentLogger::get()->info("Tabulation params: {}", tabulationParams.dump());

        jgap::Tabulate::tabulate(potential, tabulationParams, outputFilePrefix);
        jgap::CurrentLogger::get()->info("Tabulation complete");

    } catch (exception& e) {
        jgap::CurrentLogger::get()->error("Fail: {}", e.what());
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}