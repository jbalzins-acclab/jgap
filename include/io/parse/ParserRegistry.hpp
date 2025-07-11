#ifndef PARSERREGISTRY_HPP
#define PARSERREGISTRY_HPP

#include <string>
#include <map>
#include <functional>
#include <memory>
#include <nlohmann/json.hpp>

using namespace std;

#ifndef REGISTER_PARSER
#define REGISTER_PARSER(name, baseType, regType) \
struct regType##Register { \
    regType##Register() { \
        ParserRegistry<baseType>::getRegistry()[name] = [](const nlohmann::json& j){return make_shared<regType>(j);}; \
    } \
}; \
static regType##Register regType##RegisterInstance;
#endif

namespace jgap {

    template<class TBase>
    class ParserRegistry {
    public:
        //inline static map<string, function<shared_ptr<TBase>(nlohmann::json)>> registry = {};
        static auto& getRegistry() {
            static auto registry = make_shared<map<string, function<shared_ptr<TBase>(nlohmann::json)>>>();
            return *registry;
        }
        static shared_ptr<TBase> get(const nlohmann::json& params) {
            return getRegistry()[params["type"]](params);
        }
    };
}

#endif //PARSERREGISTRY_HPP
