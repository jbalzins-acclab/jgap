#include <core/fit/InRamJgapFit.hpp>
#include <utils/Utils.hpp>
#include <iostream>
#include <fstream>
#include <nlohmann/json.hpp>
#include <format>
#include <ParserRegistryAuto.hpp>
#include <Version.hpp>
#include <tbb/parallel_for.h>

#include "core/fit/Tabulate.hpp"
#include "io/convert/QuipXmlConverter.hpp"

using namespace std;

int main(int argc, char** argv) {

    jgap::CurrentLogger::get()->info(format("jGAP from QUIP xml v{}", JGAP_VERSION));

    if (argc != 2) {
        jgap::CurrentLogger::get()->error(
            "Expected 1 arguments: \n 1. <<tabulation-params>>.json"
            );
        return EXIT_FAILURE;
    }

    try {

        // ------------------------ READ PARAMS AND PREPARE -------------------------------
        jgap::CurrentLogger::get()->info(format("Converting QUIP GAP to jGAP: {}", argv[1]));

        pugi::xml_document quipDocument;
        if (!quipDocument.load_file(argv[1])) {
            jgap::CurrentLogger::get()->error(format("Cannot input quip file: {}", argv[1]));
            return EXIT_FAILURE;
        }

        jgap::CurrentLogger::get()->info("Converting");
        auto result = jgap::QuipXmlConverter::transform(quipDocument.document_element());

        string outputFilename = string(argv[1]) + ".jgap.json";
        jgap::CurrentLogger::get()->info("Converted => saving to: " + outputFilename);

        ofstream outputFile(outputFilename);
        if (!outputFile.is_open()) {
            jgap::CurrentLogger::get()->error(format("Cannot open output file {}", outputFilename));
            return EXIT_FAILURE;
        }
        outputFile << result->serialize().dump(4);
        outputFile.flush();
        outputFile.close();

    } catch (exception& e) {
        jgap::CurrentLogger::get()->error("Fail: " + string(e.what()));
        throw;
    }

    return EXIT_SUCCESS;
}