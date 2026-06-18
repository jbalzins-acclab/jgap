#include "TabGapPotentialSerialization.hpp"
#include "core/potentials/tabgap/TabGapPotential.hpp"
#include "io/tabgap/TabGapIO.hpp"
#include "serialization/SerializationNode.hpp"

namespace jgap {

    bool TabGapPotentialSerialization::serialize(const ValuePtr<Potential>& obj, SerializationNode& node) const {
        auto derived = obj.as<TabGapPotential>();
        if (!derived) {
            return false;
        }

        // Writes the tabGAP layout (incl. jgap=true and the embedded eam_files) into the node's group;
        // no separate .eam.fs files are produced.
        TabGapIO::writeToGroup(node.hdfGroup(), *derived);
        return true;
    }

    ValuePtr<Potential> TabGapPotentialSerialization::deserialize(const SerializationNode& node) const {
        // The jgap=true marker (written by TabGapIO) is what identifies a tabGAP node.
        if (node.readOptionalAttribute<bool>("jgap") != true) {
            return nullptr;
        }

        return ValuePtr<Potential>(TabGapIO::fromGroup(node.hdfGroup(), /*read_embedded_eam_fs=*/true));
    }

    REGISTER_SERIALIZATION(TabGapPotentialSerialization, Potential);
}
