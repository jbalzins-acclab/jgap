#include "TabGapPotentialSerialization.hpp"
#include "jgap/core/potentials/tabgap/TabGapPotential.hpp"
#include "jgap/io/tabgap/TabGapIO.hpp"
#include "jgap/serialization/SerializationNode.hpp"

namespace jgap {

    bool TabGapPotentialSerialization::serialize(const ValuePtr<Potential>& obj, SerializationNode& node) const {
        auto derived = obj.as<TabGapPotential>();
        if (!derived) {
            return false;
        }

        auto* h5_group = node.asH5();
        if (!h5_group) {
            JGAP_LOG_AND_THROW("SerializationNode is not an HDF5 node");
        }

        TabGapIO::writeToGroup(*h5_group, *derived);
        return true;
    }

    ValuePtr<Potential> TabGapPotentialSerialization::deserialize(const SerializationNode& node) const {
        auto* h5_group = node.asH5();
        if (!h5_group) {
            return nullptr;
        }

        TabGapPotential pot = TabGapIO::fromGroup(*h5_group, /*read_embedded_eam_fs=*/true);
        if (pot.getTwoBodyComponents().empty() && pot.getThreeBodyComponents().empty() &&
            pot.getEamComponents().empty() && pot.getIsolatedAtomEnergies().empty()) {
            return nullptr;
        }

        return ValuePtr<Potential>(std::move(pot));
    }

    REGISTER_SERIALIZATION(TabGapPotentialSerialization, Potential);
}
