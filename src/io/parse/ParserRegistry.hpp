#ifndef PARSERREGISTRY_HPP
#define PARSERREGISTRY_HPP

#include <string>
#include <map>
#include <functional>
#include <memory>

#include "utils/Utils.hpp"

/*
namespace jgap {
    template<class TBase>
    class ParserRegistry {
    public:
        static auto &getRegistry() {
            static auto registry = std::make_shared<
                std::map<std::string, std::function<std::shared_ptr<TBase>(DataNode)> >
            >();
            return *registry;
        }
    };
}
*/

#define ADD_PARSER(BaseType, RegType, type_id) \
    struct RegType##BaseType##type_id##Register { \
    RegType##BaseType##type_id##Register() { \
    jgap::ParserRegistry<BaseType>::getRegistry()[RegType::TYPE_ID] = [](const jgap::DataNode &n){ \
    return RegType::fromDataNode(n); \
    }; \
    } \
    }; \
    static RegType##BaseType##type_id##Register RegType##BaseType##type_id##RegisterInstance;

#define SETUP_PARSER(BaseType, RegType, type_id) \
    static constexpr std::string TYPE_ID = #type_id; \
    static std::shared_ptr<RegType> fromDataNode(const DataNode &params); \
    struct RegType##BaseType##type_id##Register { \
    RegType##BaseType##type_id##Register() { \
    jgap::ParserRegistry<BaseType>::getRegistry()[RegType::TYPE_ID] = [](const jgap::DataNode &n){ \
    return RegType::fromDataNode(n); \
    }; \
    } \
    }; \
    static RegType##BaseType##type_id##Register RegType##BaseType##type_id##RegisterInstance;

#define SETUP_PARSER_AND_SERIALIZATION(BaseType, RegType, type_id) \
    SETUP_PARSER(BaseType, RegType, type_id) \
    std::string getTypeId() override { return TYPE_ID; }\
    DataNode serialize() override;

#define REGISTRY_GET(BaseType, params) \
    jgap::ParserRegistry<BaseType>::getRegistry().contains(REQUIRE(params, "type")) ? \
    jgap::ParserRegistry<BaseType>::getRegistry()[REQUIRE(params, "type")](params) : \
    ([&]() -> std::shared_ptr<BaseType> { \
    JGAP_LOG_AND_THROW("type={} is not a known(registered) implementation of {}", \
    REQUIRE(params, "type").asString(), #BaseType); \
    })()

#endif
