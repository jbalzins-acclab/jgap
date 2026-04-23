#ifndef JGAP_YAMLCONVERTER_HPP
#define JGAP_YAMLCONVERTER_HPP

#include "data/DataNode.hpp"

namespace jgap {
    class YamlConverter {
    public:
        static YAML::Node toYaml(const DataNode& node);
        static DataNode fromYaml(const YAML::Node&);
    };
}

#endif