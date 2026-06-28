#ifndef JGAP_CLUSTERSYMMETRY_HPP
#define JGAP_CLUSTERSYMMETRY_HPP

namespace jgap {

    enum class ClusterSymmetry {
        IndexSensitive,
        NodeSymmetric,
        FullSymmetry
    };

    using ClusterSymmetry::IndexSensitive;
    using ClusterSymmetry::NodeSymmetric;
    using ClusterSymmetry::FullSymmetry;
}

#endif
