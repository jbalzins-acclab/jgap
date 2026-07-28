// Converts a QUIP potential .xml into a jgap .h5 (the serialized potential), written with the same
// prefix as the input. Built only when XML support (pugixml) is available; see CMakeLists.txt.
//
// Usage: jgap_convert <potential.quip.xml> [output.h5]
//   - With one argument, the output is "<prefix>.h5" where <prefix> is the input path with a trailing
//     ".xml" (and a preceding ".quip", if present) stripped, e.g. Fe.quip.xml -> Fe.h5.
//   - With two arguments, the second is used verbatim as the output path.

#include <iostream>
#include <string>

#include "jgap/io/convert/QuipXmlConverter.hpp"
#include "jgap/core/potentials/Potential.hpp"
#include "jgap/core/ValuePtr.hpp"
#include "jgap/serialization/SerializationRegistry.hpp"
#include "jgap/io/log/CurrentLogger.hpp"

using namespace jgap;

std::string defaultOutput(std::string path) {
    if (path.ends_with(".xml")) {
        path.resize(path.size() - 4);
    }
    return path + ".h5";
}

int main(int argc, char** argv) {
    CurrentLogger::initDefault({});

    if (argc < 2 || argc > 3) {
        std::cerr << "Usage: " << argv[0] << " <quip-potential.xml> [output.h5]\n";
        return 1;
    }

    const std::string input = argv[1];
    const std::string output = (argc == 3) ? std::string(argv[2]) : defaultOutput(input);

    const ValuePtr<Potential> potential = QuipXmlConverter::transform(input);

    SerializationRegistry<Potential>::serialize(potential, output);
    JGAP_LOG_INFO("Converted {} -> {}", input, output);

    return 0;
}
