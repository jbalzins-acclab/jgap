#include "jgap/serialization/SerializationNode.hpp"

#include <cassert>
#include <highfive/H5DataSet.hpp>
#include <highfive/H5File.hpp>
#include <highfive/H5Group.hpp>

#include <any>
#include <typeindex>
#include <unordered_map>

namespace jgap {
    namespace detail {
        std::any& getRegistryAny(const std::type_info& type) {
            static std::unordered_map<std::type_index, std::any> global_map;
            return global_map[std::type_index(type)];
        }
    }

    struct SerializationNodeImpl {
        std::unique_ptr<HighFive::File> file_owner;
        HighFive::Group group;

        explicit SerializationNodeImpl(std::unique_ptr<HighFive::File> file)
            : file_owner(std::move(file)), group(file_owner->getGroup("/")) {}

        explicit SerializationNodeImpl(HighFive::Group grp)
            : group(std::move(grp)) {}
    };

    SerializationNode::~SerializationNode() = default;
    SerializationNode::SerializationNode(SerializationNode&&) noexcept = default;
    SerializationNode& SerializationNode::operator=(SerializationNode&&) noexcept = default;

    SerializationNode::SerializationNode(std::unique_ptr<SerializationNodeImpl> impl)
        : impl(std::move(impl)) {}

    SerializationNode SerializationNode::create(const std::string& filename) {
        return SerializationNode(std::make_unique<SerializationNodeImpl>(
            std::make_unique<HighFive::File>(filename, HighFive::File::Overwrite)));
    }

    SerializationNode SerializationNode::open(const std::string& filename) {
        return SerializationNode(std::make_unique<SerializationNodeImpl>(
            std::make_unique<HighFive::File>(filename, HighFive::File::ReadOnly)));
    }

    HighFive::Group* SerializationNode::asH5() const {
        return impl ? &impl->group : nullptr;
    }

    std::optional<SerializationNode> SerializationNode::getGroup(const std::string& name) const {
        if (!impl->group.exist(name)) {
            return std::nullopt;
        }
        return SerializationNode(std::make_unique<SerializationNodeImpl>(impl->group.getGroup(name)));
    }

    SerializationNode SerializationNode::createGroup(const std::string& group_name) {
        if (impl->group.exist(group_name)) {
            JGAP_LOG_AND_THROW("Node {} already exists", group_name);
        }
        return SerializationNode(std::make_unique<SerializationNodeImpl>(impl->group.createGroup(group_name)));
    }

    std::vector<std::string> SerializationNode::getChildNames() const {
        return impl->group.listObjectNames();
    }

    std::vector<std::string> SerializationNode::getAttributeNames() const {
        return impl->group.listAttributeNames();
    }

    // --- Attribute Writes ---
    void SerializationNode::writeAttribute(const std::string& name, int value) {
        impl->group.createAttribute<int>(name, value);
    }
    void SerializationNode::writeAttribute(const std::string& name, long value) {
        impl->group.createAttribute<long>(name, value);
    }
    void SerializationNode::writeAttribute(const std::string& name, unsigned long value) {
        impl->group.createAttribute<unsigned long>(name, value);
    }
    void SerializationNode::writeAttribute(const std::string& name, unsigned long long value) {
        impl->group.createAttribute<unsigned long long>(name, value);
    }
    void SerializationNode::writeAttribute(const std::string& name, float value) {
        impl->group.createAttribute<float>(name, value);
    }
    void SerializationNode::writeAttribute(const std::string& name, double value) {
        impl->group.createAttribute<double>(name, value);
    }
    void SerializationNode::writeAttribute(const std::string& name, bool value) {
        impl->group.createAttribute<bool>(name, value);
    }
    void SerializationNode::writeAttribute(const std::string& name, const std::string& value) {
        impl->group.createAttribute<std::string>(name, value);
    }
    void SerializationNode::writeAttribute(const std::string& name, const char* value) {
        impl->group.createAttribute<std::string>(name, std::string(value));
    }

    // --- Attribute Reads ---
    std::optional<int> SerializationNode::readOptionalIntAttribute(const std::string& name) const {
        if (!impl->group.hasAttribute(name)) return std::nullopt;
        int val;
        impl->group.getAttribute(name).read(val);
        return val;
    }
    int SerializationNode::readIntAttribute(const std::string& name) const {
        auto opt = readOptionalIntAttribute(name);
        if (!opt) JGAP_LOG_AND_THROW("Required attribute '{}' not found.", name);
        return opt.value();
    }

    std::optional<size_t> SerializationNode::readOptionalSizeAttribute(const std::string& name) const {
        if (!impl->group.hasAttribute(name)) return std::nullopt;
        size_t val;
        impl->group.getAttribute(name).read(val);
        return val;
    }
    size_t SerializationNode::readSizeAttribute(const std::string& name) const {
        auto opt = readOptionalSizeAttribute(name);
        if (!opt) JGAP_LOG_AND_THROW("Required attribute '{}' not found.", name);
        return opt.value();
    }

    std::optional<double> SerializationNode::readOptionalDoubleAttribute(const std::string& name) const {
        if (!impl->group.hasAttribute(name)) return std::nullopt;
        double val;
        impl->group.getAttribute(name).read(val);
        return val;
    }
    double SerializationNode::readDoubleAttribute(const std::string& name) const {
        auto opt = readOptionalDoubleAttribute(name);
        if (!opt) JGAP_LOG_AND_THROW("Required attribute '{}' not found.", name);
        return opt.value();
    }

    std::optional<bool> SerializationNode::readOptionalBoolAttribute(const std::string& name) const {
        if (!impl->group.hasAttribute(name)) return std::nullopt;
        bool val;
        impl->group.getAttribute(name).read(val);
        return val;
    }
    bool SerializationNode::readBoolAttribute(const std::string& name) const {
        auto opt = readOptionalBoolAttribute(name);
        if (!opt) JGAP_LOG_AND_THROW("Required attribute '{}' not found.", name);
        return opt.value();
    }

    std::optional<std::string> SerializationNode::readOptionalStringAttribute(const std::string& name) const {
        if (!impl->group.hasAttribute(name)) return std::nullopt;
        std::string val;
        impl->group.getAttribute(name).read(val);
        return val;
    }
    std::string SerializationNode::readStringAttribute(const std::string& name) const {
        auto opt = readOptionalStringAttribute(name);
        if (!opt) JGAP_LOG_AND_THROW("Required attribute '{}' not found.", name);
        return opt.value();
    }

    // --- DataSet Writes ---
    void SerializationNode::writeDataSet(const std::string& name, const std::vector<int>& data) {
        impl->group.createDataSet(name, data);
    }
    void SerializationNode::writeDataSet(const std::string& name, const std::vector<size_t>& data) {
        impl->group.createDataSet(name, data);
    }
    void SerializationNode::writeDataSet(const std::string& name, const std::vector<Real>& data) {
        impl->group.createDataSet(name, data);
    }
    void SerializationNode::writeDataSet(const std::string& name, const std::vector<std::string>& data) {
        impl->group.createDataSet(name, data);
    }
    void SerializationNode::writeDataSet(const std::string& name, const std::string& data) {
        impl->group.createDataSet(name, data);
    }

    // --- Descriptors Writing & Reading ---
    void SerializationNode::writeDescriptors(const std::string& name, const DescriptorsView& descriptors) {
        if (descriptors.empty()) {
            HighFive::DataSpace space({0, descriptors.getDim()});
            impl->group.createDataSet<double>(name, space);
            return;
        }
        HighFive::DataSpace space({descriptors.getNumDescriptors(), descriptors.getDim()});
        impl->group.createDataSet<double>(name, space).write_raw(descriptors.rawData());
    }

    DynamicDescriptors SerializationNode::readDescriptors(const std::string& name, size_t dim) const {
        auto dataset = impl->group.getDataSet(name);
        auto dims = dataset.getSpace().getDimensions();
        assert(dims.size() == 2);
        assert(dims[1] == dim);
        std::vector<Real> flat_data(dims[0] * dims[1]);
        if (dims[0] > 0) {
            dataset.read_raw(flat_data.data());
        }
        return DynamicDescriptors{std::move(flat_data), dim};
    }

    std::optional<DynamicDescriptors> SerializationNode::readOptionalDescriptors(const std::string& name, size_t dim) const {
        if (!impl->group.exist(name)) {
            return std::nullopt;
        }
        return readDescriptors(name, dim);
    }

    // --- DataSet Reads ---
    std::optional<std::vector<Real>> SerializationNode::readOptionalRealVectorDataSet(const std::string& name) const {
        if (!impl->group.exist(name)) return std::nullopt;
        std::vector<Real> data;
        impl->group.getDataSet(name).read(data);
        return data;
    }
    std::vector<Real> SerializationNode::readRealVectorDataSet(const std::string& name) const {
        auto opt = readOptionalRealVectorDataSet(name);
        if (!opt) JGAP_LOG_AND_THROW("Required dataset '{}' not found.", name);
        return opt.value();
    }

    std::optional<std::vector<size_t>> SerializationNode::readOptionalSizeVectorDataSet(const std::string& name) const {
        if (!impl->group.exist(name)) return std::nullopt;
        std::vector<size_t> data;
        impl->group.getDataSet(name).read(data);
        return data;
    }
    std::vector<size_t> SerializationNode::readSizeVectorDataSet(const std::string& name) const {
        auto opt = readOptionalSizeVectorDataSet(name);
        if (!opt) JGAP_LOG_AND_THROW("Required dataset '{}' not found.", name);
        return opt.value();
    }

    std::optional<std::vector<std::string>> SerializationNode::readOptionalStringVectorDataSet(const std::string& name) const {
        if (!impl->group.exist(name)) return std::nullopt;
        std::vector<std::string> data;
        impl->group.getDataSet(name).read(data);
        return data;
    }
    std::vector<std::string> SerializationNode::readStringVectorDataSet(const std::string& name) const {
        auto opt = readOptionalStringVectorDataSet(name);
        if (!opt) JGAP_LOG_AND_THROW("Required dataset '{}' not found.", name);
        return opt.value();
    }

    std::optional<std::string> SerializationNode::readOptionalStringDataSet(const std::string& name) const {
        if (!impl->group.exist(name)) return std::nullopt;
        std::string data;
        impl->group.getDataSet(name).read(data);
        return data;
    }
    std::string SerializationNode::readStringDataSet(const std::string& name) const {
        auto opt = readOptionalStringDataSet(name);
        if (!opt) JGAP_LOG_AND_THROW("Required dataset '{}' not found.", name);
        return opt.value();
    }
}
