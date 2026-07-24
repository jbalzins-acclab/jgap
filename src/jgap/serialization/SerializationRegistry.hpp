#ifndef JGAP_SERIALIZATIONREGISTRY_HPP
#define JGAP_SERIALIZATIONREGISTRY_HPP

#include <memory>
#include <string>
#include <vector>

#include "Serialization.hpp"
#include "jgap/core/ValuePtr.hpp"
#include "jgap/utils/Utils.hpp"

#include <any>
#include <typeinfo>

namespace jgap {

    namespace detail {
        std::any& getRegistryAny(const std::type_info& type);
    }

    template<Cloneable TBase>
    class SerializationRegistry {
    public:
        using Registry = std::vector<std::unique_ptr<Serialization<TBase>>>;

        static Registry& getRegistry() {
            std::any& storage = detail::getRegistryAny(typeid(TBase));
            if (!storage.has_value()) {
                storage = std::make_shared<Registry>();
            }
            return *std::any_cast<std::shared_ptr<Registry>>(storage);
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

/// Registers a serializer class with the SerializationRegistry for the specified base type(s).
/// Uses `__attribute__((used))` to prevent compiler/linker dead-code elimination (especially under
/// -O3 or LTO) from stripping static registration variables in anonymous namespaces.
#define REGISTER_SERIALIZATION(Serialization, ...)                                                                    \
    namespace {                                                                                                       \
        struct CONCAT(SerializationRegister, __LINE__) {                                                              \
            CONCAT(SerializationRegister, __LINE__)() {                                                               \
                jgap::SerializationRegistry<__VA_ARGS__>::getRegistry().push_back(std::make_unique<Serialization>()); \
            }                                                                                                         \
        };                                                                                                            \
        [[maybe_unused]] __attribute__((used)) static CONCAT(SerializationRegister, __LINE__)                         \
            CONCAT(SerializationRegisterInstance, __LINE__);                                                          \
    }

#endif
