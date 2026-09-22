#ifndef JGAP_SERIALIZATION_HPP
#define JGAP_SERIALIZATION_HPP
#include "SerializationNode.hpp"
#include "jgap/core/ValuePtr.hpp"

namespace jgap {


    template<Cloneable TBase>
    class Serialization {
    public:
        virtual ~Serialization() = default;

        virtual bool serialize(const ValuePtr<TBase>& obj, SerializationNode& node) const = 0;

        virtual ValuePtr<TBase> deserialize(const SerializationNode& node) const = 0;

        void serialize(const ValuePtr<TBase>& obj, const std::string& filename) const {
            SerializationNode node = SerializationNode::create(filename);
            serialize(obj, node);
        }

        virtual ValuePtr<TBase> deserialize(const std::string& filename) const {
            const SerializationNode node = SerializationNode::open(filename);
            return deserialize(node);
        }
    };
}

#endif
