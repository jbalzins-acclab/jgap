#ifndef JGAP_SERIALIZATIONREGISTRY_HPP
#define JGAP_SERIALIZATIONREGISTRY_HPP

#include <memory>
#include <string>
#include <vector>

#include "Serialization.hpp"
#include "jgap/core/ValuePtr.hpp"
#include "jgap/utils/Utils.hpp"

#include <string_view>

namespace jgap {

    namespace detail {
        void*& getRegistryPtr(std::string_view type_name);
    }

    template<Cloneable TBase>
    class SerializationRegistry {
    public:
        using Registry = std::vector<std::unique_ptr<Serialization<TBase>>>;

        static Registry& getRegistry() {
            void*& ptr = detail::getRegistryPtr(typeid(TBase).name());
            if (!ptr) {
                ptr = new Registry();
            }
            return *static_cast<Registry*>(ptr);
        }

        static void serialize(const ValuePtr<TBase>& obj, SerializationNode& node) {
            for (const auto& serializer: getRegistry()) {
                if (serializer->serialize(obj, node)) {
                    return;
                }
            }
            JGAP_LOG_AND_THROW(
                "Could not find an appropriate serializer for {}. Check the registry.", typeName<TBase>()
            );
        }

        static ValuePtr<TBase> deserialize(const SerializationNode& node) {
            for (const auto& serializer: getRegistry()) {
                auto obj = serializer->deserialize(node);
                if (obj) {
                    return obj;
                }
            }

            JGAP_LOG_AND_THROW(
                "Could not find an appropriate deserializer for {}. Check the registry.", typeName<TBase>()
            );
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

#if defined(__GNUC__) || defined(__clang__)
#define REGISTER_SERIALIZATION(Serialization, ...)                                                            \
    __attribute__((constructor)) static void CONCAT(SerializationRegisterFunc_, __LINE__)() {                 \
        jgap::SerializationRegistry<__VA_ARGS__>::getRegistry().push_back(std::make_unique<Serialization>()); \
    }
#else
#define REGISTER_SERIALIZATION(Serialization, ...)                                                                    \
    namespace {                                                                                                       \
        struct CONCAT(SerializationRegister, __LINE__) {                                                              \
            CONCAT(SerializationRegister, __LINE__)() {                                                               \
                jgap::SerializationRegistry<__VA_ARGS__>::getRegistry().push_back(std::make_unique<Serialization>()); \
            }                                                                                                         \
        };                                                                                                            \
        [[maybe_unused]] __attribute__((used)) static CONCAT(SerializationRegister, __LINE__) CONCAT(                 \
            SerializationRegisterInstance, __LINE__                                                                   \
        );                                                                                                            \
    }
#endif

#endif
