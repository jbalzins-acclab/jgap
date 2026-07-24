#ifndef JGAP_COMPOSITEPOTENTIALSERIALIZATION_HPP
#define JGAP_COMPOSITEPOTENTIALSERIALIZATION_HPP

#include "jgap/core/potentials/Potential.hpp"
#include "jgap/serialization/Serialization.hpp"
#include "jgap/serialization/SerializationRegistry.hpp"
#include "jgap/core/ValuePtr.hpp"

namespace jgap {
    class CompositePotentialSerialization : public Serialization<Potential> {
    public:
        bool serialize(const ValuePtr<Potential>& obj, SerializationNode& node) const override;
        ValuePtr<Potential> deserialize(const SerializationNode& node) const override;
    };
}

#endif
