#ifndef JGAP_SERIALIZATIONREGISTRY_HPP
#define JGAP_SERIALIZATIONREGISTRY_HPP

#include <string>
#include <vector>
#include <memory>

#include "Serialization.hpp"
#include "core/ValuePtr.hpp"
#include "utils/Utils.hpp"

namespace jgap {

    /**
     * Global collection of {@link Serialization}s for a given base type, and the entry point for
     * (de)serializing {@code TBase} objects without knowing their concrete type.
     *
     * {@code serialize} / {@code deserialize} simply try every registered serializer in turn and use the
     * first that accepts the object (returns true) or the node (returns non-null); if none does they
     * throw. Serializers add themselves at static-initialization time via {@link REGISTER_SERIALIZATION},
     * so merely linking a serializer's translation unit makes it available here.
     *
     * @tparam TBase the Cloneable polymorphic base type being (de)serialized.
     */
    template<Cloneable TBase>
    class SerializationRegistry {
    public:
        using Registry = std::vector<std::unique_ptr<Serialization<TBase>>>;

        /**
         * The mutable list of serializers for {@code TBase}; populated by {@link REGISTER_SERIALIZATION}.
         *
         * One registry per base type. jgap is built as a shared library, so this weak function-local
         * static is coalesced to a single instance across the library/client boundary by the dynamic
         * loader; registration and lookup therefore always see the same vector.
         *
         * @return the (process-wide) registry for {@code TBase}.
         */
        static Registry& getRegistry() {
            static Registry registry;
            return registry;
        }

        /**
         * Serializes {@code obj} into {@code node} using the first registered serializer that accepts it.
         *
         * @param obj the object to serialize.
         * @param node the destination node.
         * @throws std::runtime_error if no registered serializer handles {@code obj}'s type.
         */
        static void serialize(const ValuePtr<TBase>& obj, SerializationNode& node) {
            for (const auto& serializer: getRegistry()) {
                if (serializer->serialize(obj, node)) {
                    return;
                }
            }
            JGAP_LOG_AND_THROW("Could not find an appropriate serializer for {}. Check the registry.",
                               typeName<TBase>());
        }

        /**
         * Reconstructs an object from {@code node} using the first registered serializer that recognises
         * it.
         *
         * @param node the node to read.
         * @return the reconstructed object.
         * @throws std::runtime_error if no registered serializer recognises the node.
         */
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

        /**
         * Convenience overload: serialize {@code obj} into a freshly created file.
         *
         * @param obj the object to serialize.
         * @param filename the file to (over)write.
         */
        static void serialize(const ValuePtr<TBase>& obj, const std::string& filename) {
            SerializationNode node = SerializationNode::create(filename);
            serialize(obj, node);
        }

        /**
         * Convenience overload: deserialize from an existing file.
         *
         * @param filename the file to read.
         * @return the reconstructed object.
         */
        static ValuePtr<TBase> deserialize(const std::string& filename) {
            const SerializationNode node = SerializationNode::open(filename);
            return deserialize(node);
        }
    };
}

#define CONCAT_IMPL(a, b) a##b
#define CONCAT(a, b) CONCAT_IMPL(a, b)

/**
 * Registers a {@link jgap::Serialization} so {@link jgap::SerializationRegistry} can find it. Place in a
 * .cpp at namespace scope, e.g. {@code REGISTER_SERIALIZATION(GapPotentialSerialization, Potential);}
 * (for a template serializer, alias the instantiation first: {@code using S =
 * NBodyGapComponentSerialization<...>; REGISTER_SERIALIZATION(S, GapComponent);}). The variadic
 * parameter is the base type and may itself contain commas (e.g. {@code ClusterTransformation<4, 3>}).
 *
 * The helper struct lives in an anonymous namespace, so it has internal linkage: otherwise two
 * serializers registered on the same line number in different .cpp files would generate identically
 * named namespace-scope structs (external linkage), an ODR violation that makes the linker merge them
 * and run one serializer's registration twice while dropping the other's.
 */
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