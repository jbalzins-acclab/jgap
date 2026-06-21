#ifndef JGAP_SERIALIZATIONNODE_HPP
#define JGAP_SERIALIZATIONNODE_HPP

#include <highfive/H5File.hpp>
#include <highfive/H5Group.hpp>
#include <highfive/H5DataSet.hpp>
#include <optional>
#include <memory>
#include <string>
#include <vector>
#include <array>

#include "core/atomic/Descriptor.hpp"
#include "io/log/CurrentLogger.hpp"

namespace jgap {

    class SerializationNode {
    public:
        static constexpr int FormatVersion = 1;
        static constexpr const char* FormatVersionAttribute = "jgap_format_version";

        static SerializationNode create(const std::string& filename) {
            auto file = std::make_unique<HighFive::File>(filename, HighFive::File::Overwrite);
            SerializationNode node(std::move(file));
            node.writeAttribute(FormatVersionAttribute, FormatVersion);
            return node;
        }

        static SerializationNode open(const std::string& filename) {
            auto file = std::make_unique<HighFive::File>(filename, HighFive::File::ReadOnly);
            SerializationNode node(std::move(file));
            if (const auto version = node.readOptionalAttribute<int>(FormatVersionAttribute)) {
                if (*version != FormatVersion) {
                    JGAP_LOG_WARN("'{}' was written with jgap serialization format v{}, but this build "
                                  "uses v{}; deserialization may be incorrect.", filename, *version, FormatVersion);
                }
            } else {
                JGAP_LOG_WARN("'{}' has no jgap serialization format version stamp (pre-versioning file, "
                              "or not a jgap serialization file).", filename);
            }
            return node;
        }

        template<typename T>
        void writeAttribute(const std::string& name, const T& value) {
            group.createAttribute<T>(name, value);
        }

        template<size_t N>
        void writeAttribute(const std::string& name, const char (&value)[N]) {
            group.createAttribute<std::string>(name, std::string(value));
        }

        template<typename T>
        std::optional<T> readOptionalAttribute(const std::string& name) const {
            if (!group.hasAttribute(name)) {
                return std::nullopt;
            }
            T value;
            group.getAttribute(name).read(value);
            return value;
        }

        template<typename T>
        T readAttribute(const std::string& name) const {
            auto opt = readOptionalAttribute<T>(name);
            if (!opt) {
                JGAP_LOG_AND_THROW("Required attribute '{}' not found.", name);
            }
            return opt.value();
        }

        template<typename T>
        void writeDataSet(const std::string& name, const T& data) {
            group.createDataSet(name, data);
        }

        template<typename T>
        std::optional<T> readOptionalDataSet(const std::string& name) const {
            if (!group.exist(name)) {
                return std::nullopt;
            }
            T data;
            group.getDataSet(name).read(data);
            return data;
        }

        template<typename T>
        T readDataSet(const std::string& name) const {
            auto opt = readOptionalDataSet<T>(name);
            if (!opt) {
                JGAP_LOG_AND_THROW("Required dataset '{}' not found.", name);
            }
            return opt.value();
        }

        template<size_t Dim>
        void saveSparsePoints(const std::vector<Descriptor<Dim>>& sparse_points) {
            std::vector<std::array<Real, Dim>> data;
            data.reserve(sparse_points.size());
            for (const auto& sp : sparse_points) {
                data.push_back(sp.value);
            }
            group.createDataSet("sparse_points", data);
        }

        template<size_t Dim>
        std::vector<Descriptor<Dim>> loadSparsePoints() const {
            std::vector<std::array<Real, Dim>> data;
            group.getDataSet("sparse_points").read(data);

            std::vector<Descriptor<Dim>> sparse_points;
            sparse_points.reserve(data.size());
            for (const auto& arr : data) {
                sparse_points.push_back({arr});
            }
            return sparse_points;
        }

        std::optional<SerializationNode> getGroup(const std::string& name) const {
            if (!group.exist(name)) {
                return std::nullopt;
            }
            return SerializationNode(group.getGroup(name));
        }

        SerializationNode createGroup(const std::string& group_name) {
            if (group.exist(group_name)) {
                JGAP_LOG_AND_THROW("Node {} already exists", group_name);
            }
            return SerializationNode(group.createGroup(group_name));
        }

        std::vector<std::string> getChildNames() const {
            return group.listObjectNames();
        }

        std::vector<std::string> getAttributeNames() const {
            return group.listAttributeNames();
        }

        HighFive::Group& hdfGroup() { return group; }

        const HighFive::Group& hdfGroup() const { return group; }

    private:
        std::unique_ptr<HighFive::File> file_owner;
        HighFive::Group group;

        explicit SerializationNode(std::unique_ptr<HighFive::File> file)
            : file_owner(std::move(file)), group(file_owner->getGroup("/")) {}

        explicit SerializationNode(HighFive::Group group) : group(std::move(group)) {}
    };
}

#endif