#ifndef JGAP_TABGAPPOTENTIALSERIALIZATION_HPP
#define JGAP_TABGAPPOTENTIALSERIALIZATION_HPP

#include "core/potentials/Potential.hpp"
#include "serialization/Serialization.hpp"
#include "serialization/SerializationRegistry.hpp"
#include "core/ValuePtr.hpp"

namespace jgap {
    // Serializes a TabGapPotential by delegating to TabGapIO, writing the externally-defined tabGAP
    // layout straight into the node's HDF5 group (with the EAM part embedded as .eam.fs datasets and a
    // jgap=true marker). Unlike TabGapIO::write it does not emit separate .eam.fs files, and on
    // deserialization it always takes the EAM part from the embedded datasets.
    class TabGapPotentialSerialization : public Serialization<Potential> {
    public:
        bool serialize(const ValuePtr<Potential>& obj, SerializationNode& node) const override;
        ValuePtr<Potential> deserialize(const SerializationNode& node) const override;
    };
}

#endif
