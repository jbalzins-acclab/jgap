#include <utils/Utils.hpp>
#include <fstream>
#include <nlohmann/json.hpp>
#include <ParserRegistryAuto.hpp>

#include "io/convert/QuipXmlConverter.hpp"

#ifndef JGAP_CONVERTQUIPXML_HPP
#define JGAP_CONVERTQUIPXML_HPP

using namespace std;

namespace jgap {
    shared_ptr<Potential> convert(string xmlFileName);
}

#endif
