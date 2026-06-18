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
    /**
     * A node in a serialized tree, backed by an HDF5 group: scalar/array values are stored as attributes
     * or datasets, and nested objects as child groups. This is the storage abstraction that
     * {@link Serialization} works against; a node is either the root of a file (see {@link #create} /
     * {@link #open}) or a child group of another node ({@link #createGroup} / {@link #getGroup}).
     *
     * Conventions: small fixed values go in attributes ({@link #writeAttribute}/{@link #readAttribute}),
     * bulk arrays in datasets ({@link #writeDataSet}/{@link #readDataSet}); the {@code readOptional*}
     * variants return {@code std::nullopt} when absent, which serializers use to detect "not my node".
     */
    class SerializationNode {
    public:
        /**
         * Opens a file for writing (overwriting any existing file) and returns its root node.
         * @param filename the file to create/overwrite.
         * @return the root node of the new file.
         */
        static SerializationNode create(const std::string& filename) {
            auto file = std::make_unique<HighFive::File>(filename, HighFive::File::Overwrite);
            return SerializationNode(std::move(file));
        }

        /**
         * Opens an existing file read-only and returns its root node.
         * @param filename the file to open.
         * @return the root node of the file.
         */
        static SerializationNode open(const std::string& filename) {
            auto file = std::make_unique<HighFive::File>(filename, HighFive::File::ReadOnly);
            return SerializationNode(std::move(file));
        }

        /**
         * Writes a scalar/array value as an attribute on this node.
         * @tparam T value type (e.g. Real, int, std::string, std::vector<std::string>).
         * @param name attribute name.
         * @param value value to store.
         */
        template<typename T>
        void writeAttribute(const std::string& name, const T& value) {
            group.createAttribute<T>(name, value);
        }

        /**
         * Overload for string literals: they decay to {@code const char[N]}, which HighFive would
         * otherwise store as an N-element char array rather than a scalar string (so
         * {@code readAttribute<std::string>} would fail). Routes them through {@code std::string} so they
         * round-trip.
         * @param name attribute name.
         * @param value string literal to store.
         */
        template<size_t N>
        void writeAttribute(const std::string& name, const char (&value)[N]) {
            group.createAttribute<std::string>(name, std::string(value));
        }

        /**
         * Reads an attribute if present.
         * @tparam T value type to read.
         * @param name attribute name.
         * @return the value, or {@code std::nullopt} if the attribute does not exist.
         */
        template<typename T>
        std::optional<T> readOptionalAttribute(const std::string& name) const {
            if (!group.hasAttribute(name)) {
                return std::nullopt;
            }
            T value;
            group.getAttribute(name).read(value);
            return value;
        }

        /**
         * Reads a required attribute.
         * @tparam T value type to read.
         * @param name attribute name.
         * @return the value.
         * @throws std::runtime_error if the attribute is missing.
         */
        template<typename T>
        T readAttribute(const std::string& name) const {
            auto opt = readOptionalAttribute<T>(name);
            if (!opt) {
                JGAP_LOG_AND_THROW("Required attribute '{}' not found.", name);
            }
            return opt.value();
        }

        /**
         * Writes bulk data as a dataset on this node.
         * @tparam T data type (e.g. std::vector<Real>, std::array<Real, N>).
         * @param name dataset name.
         * @param data data to store.
         */
        template<typename T>
        void writeDataSet(const std::string& name, const T& data) {
            group.createDataSet(name, data);
        }

        /**
         * Reads a dataset if present.
         * @tparam T data type to read into.
         * @param name dataset name.
         * @return the data, or {@code std::nullopt} if the dataset does not exist.
         */
        template<typename T>
        std::optional<T> readOptionalDataSet(const std::string& name) const {
            if (!group.exist(name)) {
                return std::nullopt;
            }
            T data;
            group.getDataSet(name).read(data);
            return data;
        }

        /**
         * Reads a required dataset.
         * @tparam T data type to read into.
         * @param name dataset name.
         * @return the data.
         * @throws std::runtime_error if the dataset is missing.
         */
        template<typename T>
        T readDataSet(const std::string& name) const {
            auto opt = readOptionalDataSet<T>(name);
            if (!opt) {
                JGAP_LOG_AND_THROW("Required dataset '{}' not found.", name);
            }
            return opt.value();
        }

        /**
         * Stores a list of descriptors as a "sparse_points" dataset (one row per descriptor value).
         * @tparam Dim descriptor dimension.
         * @param sparse_points the descriptors to store.
         */
        template<size_t Dim>
        void saveSparsePoints(const std::vector<Descriptor<Dim>>& sparse_points) {
            std::vector<std::array<Real, Dim>> data;
            data.reserve(sparse_points.size());
            for (const auto& sp : sparse_points) {
                data.push_back(sp.value);
            }
            group.createDataSet("sparse_points", data);
        }

        /**
         * Loads descriptors previously written by {@link #saveSparsePoints}.
         * @tparam Dim descriptor dimension.
         * @return the stored descriptors.
         */
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

        /**
         * Returns a child node if the named child group exists.
         * @param name child group name.
         * @return the child node, or {@code std::nullopt} if it does not exist.
         */
        std::optional<SerializationNode> getGroup(const std::string& name) const {
            if (!group.exist(name)) {
                return std::nullopt;
            }
            return SerializationNode(group.getGroup(name));
        }

        /**
         * Creates and returns a new child node.
         * @param group_name child group name.
         * @return the new child node.
         * @throws std::runtime_error if a child of that name already exists.
         */
        SerializationNode createGroup(const std::string& group_name) {
            if (group.exist(group_name)) {
                JGAP_LOG_AND_THROW("Node {} already exists", group_name);
            }
            return SerializationNode(group.createGroup(group_name));
        }

        /** @return the names of this node's child groups/datasets. */
        std::vector<std::string> getChildNames() const {
            return group.listObjectNames();
        }

        /** @return the names of this node's attributes. */
        std::vector<std::string> getAttributeNames() const {
            return group.listAttributeNames();
        }

        /**
         * Direct access to the underlying HDF5 group, for serializers that must follow an externally
         * defined on-disk layout (e.g. the tabGAP format handled by TabGapIO) rather than this node's
         * generic attribute/dataset helpers.
         * @return the backing HDF5 group.
         */
        HighFive::Group& hdfGroup() { return group; }

        /** @return the backing HDF5 group (const). */
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