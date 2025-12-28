#include "io/convert/YamlConverter.hpp"

namespace jgap {
    YAML::Node YamlConverter::toYaml(const DataNode& node) {
        switch (node.type) {
            case DataNode::Type::BOOL: return YAML::Node(std::get<bool>(node.value));
            case DataNode::Type::Int: return YAML::Node(std::get<int64_t>(node.value));
            case DataNode::Type::Double: return YAML::Node(std::get<double>(node.value));
            case DataNode::Type::String: return YAML::Node(std::get<std::string>(node.value));
            case DataNode::Type::Array: {
                YAML::Node seq(YAML::NodeType::Sequence);
                for (const auto& x : std::get<std::vector<DataNode>>(node.value)) {
                    seq.push_back(toYaml(x));
                }
                return seq;
            }
            case DataNode::Type::Object: {
                YAML::Node map(YAML::NodeType::Map);
                for (const auto& [k,v] : std::get<std::map<std::string, DataNode>>(node.value)) {
                    map[k] = toYaml(v);
                }
                return map;
            }
        }
        return {};
    }

    DataNode YamlConverter::fromYaml(const YAML::Node &yaml) {
        if (!yaml) return DataNode::object();
        if (yaml.IsScalar()) {
            auto s = yaml.as<std::string>();
            // Try to detect simple types
            if (s == "true") return true;
            if (s == "false") return false;
            try { return std::stoll(s); } catch (...) {}
            try { return std::stod(s); } catch (...) {}
            return s;
        }
        if (yaml.IsMap()) {
            DataNode obj = DataNode::object();
            auto& m = std::get<std::map<std::string, DataNode>>(obj.value);
            for (auto it : yaml) {
                m[it.first.as<std::string>()] = fromYaml(it.second);
            }
            return obj;
        }
        if (yaml.IsSequence()) {
            DataNode arr = DataNode::array();
            auto& v = std::get<std::vector<DataNode>>(arr.value);
            for (auto item : yaml) v.push_back(fromYaml(item));
            return arr;
        }
        return DataNode::object();
    }
}
