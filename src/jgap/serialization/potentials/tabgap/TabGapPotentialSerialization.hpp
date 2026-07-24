#ifndef JGAP_TABGAPPOTENTIALSERIALIZATION_HPP
#define JGAP_TABGAPPOTENTIALSERIALIZATION_HPP

#include "jgap/core/potentials/Potential.hpp"
#include "jgap/serialization/Serialization.hpp"
#include "jgap/serialization/SerializationRegistry.hpp"
#include "jgap/core/ValuePtr.hpp"

namespace jgap {
    // Serializes a TabGapPotential by delegating to TabGapIO, writing the externally-defined tabGAP
    // layout straight into the node (with the EAM part embedded as .eam.fs datasets).
    // Unlike TabGapIO::write it does not emit separate .eam.fs files, and on
    // deserialization it takes the EAM part from embedded datasets if present.
    class TabGapPotentialSerialization : public Serialization<Potential> {
    public:
        bool serialize(const ValuePtr<Potential>& obj, SerializationNode& node) const override;
        ValuePtr<Potential> deserialize(const SerializationNode& node) const override;
    };
}

#endif
