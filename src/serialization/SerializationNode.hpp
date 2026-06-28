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
        static SerializationNode create(const std::string& filename) {
            return SerializationNode(std::make_unique<HighFive::File>(filename, HighFive::File::Overwrite));
        }

        static SerializationNode open(const std::string& filename) {
            return SerializationNode(std::make_unique<HighFive::File>(filename, HighFive::File::ReadOnly));
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