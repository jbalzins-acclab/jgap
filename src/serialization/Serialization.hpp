#ifndef JGAP_SERIALIZATION_HPP
#define JGAP_SERIALIZATION_HPP
#include "SerializationNode.hpp"
#include "core/ValuePtr.hpp"

namespace jgap {

    /**
     * Strategy for (de)serializing one concrete subtype of the polymorphic base {@code TBase} to/from a
     * {@link SerializationNode} (an HDF5 group). There is one Serialization per concrete type; they are
     * collected in {@link SerializationRegistry} which tries each in turn, so an individual Serialization
     * must recognise "its own" objects/nodes and politely decline the rest.
     *
     * <h3>Adding a serializer for a new type</h3>
     *
     * <b>1. The common case</b> &mdash; a concrete subtype of a Cloneable polymorphic base. Derive from
     * {@code Serialization<TBase>} and:
     * <ul>
     *   <li>in {@code serialize}, dynamic-cast the object to your concrete type with
     *       {@code obj.as<Concrete>()}; if it is not yours return {@code false} so the registry tries the
     *       next serializer. Otherwise write a discriminator (by convention a {@code "name"} string
     *       attribute) plus the object's data, and return {@code true}.</li>
     *   <li>in {@code deserialize}, check that discriminator with
     *       {@code node.readOptionalAttribute<std::string>("name")}; if it does not match return
     *       {@code nullptr} (again so the registry moves on), otherwise rebuild and return the object.</li>
     *   <li>register the serializer with {@code REGISTER_SERIALIZATION(MySerialization, TBase)} in the
     *       .cpp (see {@link SerializationRegistry}).</li>
     * </ul>
     * See IsolatedAtomPotentialSerialization for a minimal example. To serialize a type that owns other
     * polymorphic objects, recurse through {@code SerializationRegistry<ChildBase>} on a child node
     * created with {@code node.createGroup(...)} / read back with {@code node.getGroup(...)}; see
     * GapPotentialSerialization (external potential + components) or CompositePotentialSerialization.
     *
     * <b>2. Template types</b> &mdash; when the serialized type is a class template instantiated for
     * several parameter sets, the pattern differs: make the serializer itself a template and, in the
     * .cpp, explicitly instantiate AND register one specialization per instantiation in use. See
     * NBodyGapComponentSerialization (and TransformationAggregatorImplSerialization).
     *
     * @tparam TBase the Cloneable polymorphic base type whose subtypes this serializes; objects are
     *               handled as {@link ValuePtr}&lt;TBase&gt;.
     */
    template<Cloneable TBase>
    class Serialization {
    public:
        virtual ~Serialization() = default;

        /**
         * Writes {@code obj} into {@code node} if it is of this serializer's concrete type.
         *
         * @param obj the object to serialize.
         * @param node the destination node (HDF5 group).
         * @return true if this serializer handled {@code obj}; false if {@code obj} is a different type
         *         and the caller should try another serializer.
         */
        virtual bool serialize(const ValuePtr<TBase>& obj, SerializationNode& node) const = 0;

        /**
         * Rebuilds an object from {@code node} if the node holds this serializer's concrete type.
         *
         * @param node the node (HDF5 group) to read.
         * @return the reconstructed object, or {@code nullptr} if {@code node} is not this type (so the
         *         caller should try another serializer).
         */
        virtual ValuePtr<TBase> deserialize(const SerializationNode& node) const = 0;

        /**
         * Convenience overload: serialize {@code obj} into a freshly created file.
         *
         * @param obj the object to serialize.
         * @param filename the file to (over)write.
         */
        void serialize(const ValuePtr<TBase>& obj, const std::string& filename) const {
            SerializationNode node = SerializationNode::create(filename);
            serialize(obj, node);
        }

        /**
         * Convenience overload: deserialize from an existing file.
         *
         * @param filename the file to read.
         * @return the reconstructed object.
         */
        virtual ValuePtr<TBase> deserialize(const std::string& filename) const {
            const SerializationNode node = SerializationNode::open(filename);
            return deserialize(node);
        }
    };
}

#endif