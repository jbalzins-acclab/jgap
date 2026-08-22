#ifndef JGAP_SCREENEDCOULOMBPOTENTIALSERIALIZATION_HPP
#define JGAP_SCREENEDCOULOMBPOTENTIALSERIALIZATION_HPP

#include "jgap/core/ValuePtr.hpp"
#include "jgap/core/potentials/Potential.hpp"
#include "jgap/serialization/Serialization.hpp"
#include "jgap/serialization/SerializationRegistry.hpp"

namespace jgap {
    class ScreenedCoulombPotentialSerialization : public Serialization<Potential> {
    public:
        bool serialize(const ValuePtr<Potential>& obj, SerializationNode& node) const override;
        ValuePtr<Potential> deserialize(const SerializationNode& node) const override;
    };
}

#endif
