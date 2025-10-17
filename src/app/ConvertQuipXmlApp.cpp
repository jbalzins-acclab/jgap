#include <utils/Utils.hpp>
#include <fstream>
#include <nlohmann/json.hpp>
#include <ParserRegistryAuto.hpp>

#include "io/convert/QuipXmlConverter.hpp"

using namespace std;

namespace jgap {

    shared_ptr<Potential> convert(string quipFileName) {
        // ------------------------ READ PARAMS AND PREPARE -------------------------------
        CurrentLogger::get()->info("Converting QUIP {} to jGAP", quipFileName);

        pugi::xml_document quipDocument;
        if (!quipDocument.load_file(quipFileName.c_str())) {
            CurrentLogger::get()->logAndThrow("Cannot input quip file: {}", quipFileName);
        }

        auto result = QuipXmlConverter::transform(quipDocument.document_element());

        string outputFilename = quipFileName + ".jgap.json";
        if (quipFileName.contains(".xml")) {
            outputFilename = quipFileName.replace(quipFileName.rfind(".xml"), 4, ".jgap.json");
        }
        CurrentLogger::get()->info("Converted => saving to: {}", outputFilename);

        ofstream outputFile(outputFilename);
        if (!outputFile.is_open()) {
            CurrentLogger::get()->logAndThrow("Could not open output file {}", outputFilename);
        }
        outputFile << result->serialize().dump(4);
        outputFile.flush();
        outputFile.close();

        return result;
    }
}