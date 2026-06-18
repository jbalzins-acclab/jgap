#ifndef JGAP_SERIALIZATIONREGISTRY_HPP
#define JGAP_SERIALIZATIONREGISTRY_HPP

#include <string>
#include <map>
#include <vector>
#include <functional>
#include <memory>
#include <typeindex>

#include "Serialization.hpp"
#include "core/ValuePtr.hpp"
#include "utils/Utils.hpp"

namespace jgap {

    namespace detail {
        // Single, externally-linked storage for all registries, keyed by base type. Defined once in
        // SerializationRegistry.cpp. Routing every registry through this one map avoids relying on the
        // coalescing of per-translation-unit weak statics (which is unreliable across the static-library
        // boundary on this toolchain), so registration and lookup always see the same vector.
        std::map<std::type_index, std::shared_ptr<void>>& serializationRegistryStorage();
    }

    template<Cloneable TBase>
    class SerializationRegistry {
    public:
        using Registry = std::vector<std::unique_ptr<Serialization<TBase>>>;

        static Registry& getRegistry() {
            auto& storage = detail::serializationRegistryStorage();
            auto key = std::type_index(typeid(TBase));
            auto it = storage.find(key);
            if (it == storage.end()) {
                auto registry = std::make_shared<Registry>();
                storage.emplace(key, registry);
                return *registry;
            }
            return *std::static_pointer_cast<Registry>(it->second);
        }

        static void serialize(const ValuePtr<TBase>& obj, SerializationNode& node) {
            for (const auto& serializer: getRegistry()) {
                if (serializer->serialize(obj, node)) {
                    return;
                }
            }
            JGAP_LOG_AND_THROW("Could not find an appropriate serializer for {}. Check the registry.",
                               typeName<TBase>());
        }

        static ValuePtr<TBase> deserialize(const SerializationNode& node) {
            for (const auto& serializer: getRegistry()) {
                auto obj = serializer->deserialize(node);
                if (obj) {
                    return obj;
                }
            }

            JGAP_LOG_AND_THROW("Could not find an appropriate deserializer for {}. Check the registry.",
                               typeName<TBase>());
        }

        static void serialize(const ValuePtr<TBase>& obj, const std::string& filename) {
            SerializationNode node = SerializationNode::create(filename);
            serialize(obj, node);
        }

        static ValuePtr<TBase> deserialize(const std::string& filename) {
            const SerializationNode node = SerializationNode::open(filename);
            return deserialize(node);
        }
    };
}

#define CONCAT_IMPL(a, b) a##b
#define CONCAT(a, b) CONCAT_IMPL(a, b)

// The helper struct lives in an anonymous namespace so it has internal linkage: otherwise two
// serializers registered on the same line number in different .cpp files would generate identically
// named namespace-scope structs (external linkage), an ODR violation that makes the linker merge them
// and run one serializer's registration twice while dropping the other's.
#define REGISTER_SERIALIZATION(Serialization, ...) \
    namespace { \
        struct CONCAT(SerializationRegister, __LINE__) { \
            CONCAT(SerializationRegister, __LINE__)() { \
                jgap::SerializationRegistry<__VA_ARGS__>::getRegistry().push_back(std::make_unique<Serialization>()); \
            } \
        }; \
        static CONCAT(SerializationRegister, __LINE__) CONCAT(SerializationRegisterInstance, __LINE__); \
    }

#endif