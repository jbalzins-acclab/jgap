#ifndef JGAP_ZBLPOTENTIALSERIALIZATION_HPP
#define JGAP_ZBLPOTENTIALSERIALIZATION_HPP

#include "core/potentials/Potential.hpp"
#include "serialization/Serialization.hpp"
#include "serialization/SerializationRegistry.hpp"
#include "core/ValuePtr.hpp"

namespace jgap {
    class ZblPotentialSerialization : public Serialization<Potential> {
    public:
        bool serialize(const ValuePtr<Potential>& obj, SerializationNode& node) const override;
        ValuePtr<Potential> deserialize(const SerializationNode& node) const override;
    };
}

#endif
