#include "io/params/JsonParamReader.hpp"

#include <fstream>
#include <exception>
#include <memory>

#include <nlohmann/json.hpp>

using json = nlohmann::json;

using namespace std;

namespace jgap {
    /*
    shared_ptr<ParamNode> JsonParamReader::read(const string_view fileName) {
        ifstream inFile{string(fileName)};
        if (!inFile) {
            throw runtime_error("Failed to open file: " + string(fileName));
        }

        const json rootNode = json::parse(inFile);

        return nullptr;//make_shared<NLJsonWrapperNode>(make_shared<json>(rootNode));
    }*/

}
