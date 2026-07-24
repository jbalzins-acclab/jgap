#ifndef JGAP_GAPPOTENTIALSERIALIZATION_HPP
#define JGAP_GAPPOTENTIALSERIALIZATION_HPP

#include "jgap/core/ValuePtr.hpp"
#include "jgap/core/potentials/Potential.hpp"
#include "jgap/serialization/Serialization.hpp"

namespace jgap {
    class GapPotentialSerialization : public Serialization<Potential> {
    public:
        bool serialize(const ValuePtr<Potential>& obj, SerializationNode& node) const override;
        ValuePtr<Potential> deserialize(const SerializationNode& node) const override;
    };
}

#endif
