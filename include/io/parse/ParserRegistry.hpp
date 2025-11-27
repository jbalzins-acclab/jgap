#ifndef PARSERREGISTRY_HPP
#define PARSERREGISTRY_HPP

#include <string>
#include <map>
#include <functional>
#include <memory>
#include <nlohmann/json.hpp>


#define REGISTER_PARSER(baseType, regType) \
struct regType##Register { \
    regType##Register() { \
        jgap::ParserRegistry<baseType>::getRegistry()[regType::TYPE] = [](const nlohmann::json& j){ \
            return regType::fromJson(j); \
        }; \
    } \
}; \
static regType##Register regType##RegisterInstance;

namespace jgap {
    template<class TBase>
    class ParserRegistry {
    public:
        static auto& getRegistry() {
            static auto registry = std::make_shared<
                std::map<std::string, std::function<std::shared_ptr<TBase>(nlohmann::json)>>
            >();
            return *registry;
        }
    };
}

#define REGISTRY_GET(baseType, params) \
    jgap::ParserRegistry<baseType>::getRegistry().contains(require(params, "type")) ? \
    jgap::ParserRegistry<baseType>::getRegistry()[require(params, "type")](params) : \
    ([&]() -> std::shared_ptr<baseType> { \
    JGAP_LOG_AND_THROW("type={} is not a known(registered) implementation of {}", \
    require(params, "type").get<std::string>(), #baseType); \
    })()

#endif
