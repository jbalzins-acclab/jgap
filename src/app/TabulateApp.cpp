#include <core/fit/QRGapFit.hpp>
#include <utils/Utils.hpp>
#include <iostream>
#include <fstream>
#include <nlohmann/json.hpp>
#include <format>
#include <ParserRegistryAuto.hpp>
#include <Version.hpp>
#include <tbb/parallel_for.h>

#include "io/convert/QuipXmlConverter.hpp"


int main(int argc, char** argv) {

    jgap::JGAP_LOG_INFO("jGAP tabulate v{}", JGAP_VERSION);

    if (argc != 2) {
        jgap::JGAP_LOG_ERROR(
            "Expected 1 arguments: \n 1. <<tabulation-params>>.json"
            );
        return EXIT_FAILURE;
    }

    try {

        // ------------------------ READ PARAMS AND PREPARE -------------------------------
        jgap::JGAP_LOG_INFO("Tabulation as specified in: {}", argv[1]);

        std::string paramFileName = argv[1];
        ifstream paramFile(paramFileName);
        if (!paramFile.is_open()) {
            jgap::JGAP_LOG_ERROR("Cannot open tabulation-param file {}", paramFileName);
            return EXIT_FAILURE;
        }

        DataNode tabulationParams;
        paramFile >> tabulationParams;

        const std::string outputFilePrefix = tabulationParams["output_file_prefix"];
        const std::string potentialFileName = tabulationParams["potential_file"];

        std::shared_ptr<jgap::Potential> potential;
        if (potentialFileName.ends_with(".xml")) {

            pugi::xml_document quipDocument;
            if (!quipDocument.load_file(potentialFileName.c_str())) {
                jgap::JGAP_LOG_ERROR("Cannot input quip file: {}",  potentialFileName);
                return EXIT_FAILURE;
            }

            jgap::JGAP_LOG_INFO("Converting from QUIP xml");
            potential = jgap::QuipXmlConverter::transform(quipDocument.document_element());

        } else {
            ifstream potentialParamFile(potentialFileName);
            if (!potentialParamFile.is_open()) {
                jgap::JGAP_LOG_ERROR("Cannot open pot-param file {}", potentialFileName);
                return EXIT_FAILURE;
            }
            DataNode potParams;
            potentialParamFile >> potParams;
            potential = jgap::ParserRegistry<jgap::Potential>::get(potParams);
        }

        jgap::JGAP_LOG_INFO("Tabulating potential: {}", potential->serialize().dump());
        jgap::JGAP_LOG_INFO("Tabulation params: {}", tabulationParams.dump());

        jgap::Tabulate::tabulate(potential, tabulationParams, outputFilePrefix);
        jgap::JGAP_LOG_INFO("Tabulation complete");

    } catch (exception& e) {
        jgap::JGAP_LOG_ERROR("Fail: {}", e.what());
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}