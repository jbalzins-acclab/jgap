// Converts a QUIP potential .xml into a jgap .h5 (the serialized potential), written with the same
// prefix as the input. Built only when XML support (pugixml) is available; see CMakeLists.txt.
//
// Usage: jgap_convert <potential.quip.xml> [output.h5]
//   - With one argument, the output is "<prefix>.h5" where <prefix> is the input path with a trailing
//     ".xml" (and a preceding ".quip", if present) stripped, e.g. Fe.quip.xml -> Fe.h5.
//   - With two arguments, the second is used verbatim as the output path.

#include <iostream>
#include <string>

#include <pugixml.hpp>

#include "io/convert/QuipXmlConverter.hpp"
#include "core/potentials/Potential.hpp"
#include "core/ValuePtr.hpp"
#include "serialization/SerializationRegistry.hpp"
#include "io/log/CurrentLogger.hpp"

using namespace jgap;

namespace {
    std::string defaultOutput(std::string path) {
        if (path.ends_with(".xml")) {
            path.resize(path.size() - 4);
        }
        if (path.ends_with(".quip")) {
            path.resize(path.size() - 5);
        }
        return path + ".h5";
    }
}

int main(int argc, char** argv) {
    CurrentLogger::initDefault({});

    if (argc < 2 || argc > 3) {
        std::cerr << "Usage: " << argv[0] << " <potential.quip.xml> [output.h5]\n";
        return 1;
    }

    const std::string input = argv[1];
    const std::string output = (argc == 3) ? std::string(argv[2]) : defaultOutput(input);

    pugi::xml_document doc;
    const pugi::xml_parse_result parsed = doc.load_file(input.c_str());
    if (!parsed) {
        std::cerr << "Failed to parse '" << input << "': " << parsed.description() << "\n";
        return 1;
    }

    const ValuePtr<Potential> potential = QuipXmlConverter::transform(doc.document_element());

    SerializationRegistry<Potential>::serialize(potential, output);
    JGAP_LOG_INFO("Converted {} -> {}", input, output);

    return 0;
}
