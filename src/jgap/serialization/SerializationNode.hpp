#ifndef JGAP_SERIALIZATIONNODE_HPP
#define JGAP_SERIALIZATIONNODE_HPP

#include <array>
#include <memory>
#include <optional>
#include <string>
#include <vector>

#include "jgap/core/Real.hpp"
#include "jgap/core/atomic/descriptor/Descriptor.hpp"
#include "jgap/core/atomic/descriptor/DescriptorsView.hpp"
#include "jgap/core/atomic/descriptor/DynamicDescriptors.hpp"
#include "jgap/io/log/CurrentLogger.hpp"

namespace HighFive {
    class Group;
}

namespace jgap {

    struct SerializationNodeImpl;

    class SerializationNode {
    public:
        ~SerializationNode();

        SerializationNode(SerializationNode&&) noexcept;
        SerializationNode& operator=(SerializationNode&&) noexcept;

        SerializationNode(const SerializationNode&) = delete;
        SerializationNode& operator=(const SerializationNode&) = delete;

        static SerializationNode create(const std::string& filename);
        static SerializationNode open(const std::string& filename);

        // --- Attribute Writes ---
        void writeAttribute(const std::string& name, int value);
        void writeAttribute(const std::string& name, long value);
        void writeAttribute(const std::string& name, unsigned long value);
        void writeAttribute(const std::string& name, unsigned long long value);
        void writeAttribute(const std::string& name, float value);
        void writeAttribute(const std::string& name, double value);
        void writeAttribute(const std::string& name, bool value);
        void writeAttribute(const std::string& name, const std::string& value);
        void writeAttribute(const std::string& name, const char* value);

        // --- Attribute Reads ---
        std::optional<int> readOptionalIntAttribute(const std::string& name) const;
        int readIntAttribute(const std::string& name) const;

        std::optional<size_t> readOptionalSizeAttribute(const std::string& name) const;
        size_t readSizeAttribute(const std::string& name) const;

        std::optional<double> readOptionalDoubleAttribute(const std::string& name) const;
        double readDoubleAttribute(const std::string& name) const;

        std::optional<bool> readOptionalBoolAttribute(const std::string& name) const;
        bool readBoolAttribute(const std::string& name) const;

        std::optional<std::string> readOptionalStringAttribute(const std::string& name) const;
        std::string readStringAttribute(const std::string& name) const;

        // --- DataSet Writes ---
        void writeDataSet(const std::string& name, const std::vector<int>& data);
        void writeDataSet(const std::string& name, const std::vector<size_t>& data);
        void writeDataSet(const std::string& name, const std::vector<double>& data);
        void writeDataSet(const std::string& name, const std::vector<std::string>& data);
        void writeDataSet(const std::string& name, const std::string& data);

        template<typename T, size_t N>
        void writeDataSet(const std::string& name, const std::array<T, N>& data) {
            writeDataSet(name, std::vector<T>(data.begin(), data.end()));
        }

        // --- Descriptors Writing & Reading ---
        void writeDescriptors(const std::string& name, const DescriptorsView& descriptors);

        template<size_t Dim>
        void writeDescriptors(const std::string& name, const std::vector<Descriptor<Dim>>& descriptors) {
            writeDescriptors(name, DescriptorsView(descriptors));
        }

        template<size_t Dim>
        void writeDataSet(const std::string& name, const std::vector<Descriptor<Dim>>& descriptors) {
            writeDescriptors(name, DescriptorsView(descriptors));
        }

        DynamicDescriptors readDescriptors(const std::string& name, size_t dim) const;
        std::optional<DynamicDescriptors> readOptionalDescriptors(const std::string& name, size_t dim) const;

        template<size_t Dim>
        std::vector<Descriptor<Dim>> readDescriptors(const std::string& name) const {
            return readDescriptors(name, Dim).template releaseAsVector<Dim>();
        }

        template<size_t Dim>
        std::optional<std::vector<Descriptor<Dim>>> readOptionalDescriptors(const std::string& name) const {
            auto opt = readOptionalDescriptors(name, Dim);
            if (!opt) return std::nullopt;
            return std::move(opt.value()).template releaseAsVector<Dim>();
        }

        // --- DataSet Reads ---
        std::optional<std::vector<double>> readOptionalDoubleVectorDataSet(const std::string& name) const;
        std::vector<double> readDoubleVectorDataSet(const std::string& name) const;

        template<size_t N>
        std::optional<std::array<Real, N>> readOptionalDoubleArrayDataSet(const std::string& name) const {
            auto opt = readOptionalDoubleVectorDataSet(name);
            if (!opt) return std::nullopt;
            const auto& vec = opt.value();
            if (vec.size() != N) {
                JGAP_LOG_AND_THROW("Dataset '{}' expected size {}, got {}", name, N, vec.size());
            }
            std::array<Real, N> arr;
            std::copy(vec.begin(), vec.end(), arr.begin());
            return arr;
        }

        template<size_t N>
        std::array<Real, N> readDoubleArrayDataSet(const std::string& name) const {
            auto opt = readOptionalDoubleArrayDataSet<N>(name);
            if (!opt) {
                JGAP_LOG_AND_THROW("Required dataset '{}' not found.", name);
            }
            return opt.value();
        }

        std::optional<std::vector<size_t>> readOptionalSizeVectorDataSet(const std::string& name) const;
        std::vector<size_t> readSizeVectorDataSet(const std::string& name) const;

        std::optional<std::vector<std::string>> readOptionalStringVectorDataSet(const std::string& name) const;
        std::vector<std::string> readStringVectorDataSet(const std::string& name) const;

        std::optional<std::string> readOptionalStringDataSet(const std::string& name) const;
        std::string readStringDataSet(const std::string& name) const;

        std::optional<SerializationNode> getGroup(const std::string& name) const;
        SerializationNode createGroup(const std::string& group_name);

        std::vector<std::string> getChildNames() const;
        std::vector<std::string> getAttributeNames() const;

        HighFive::Group* asH5() const;

    private:
        explicit SerializationNode(std::unique_ptr<SerializationNodeImpl> impl);

        std::unique_ptr<SerializationNodeImpl> impl;
    };
}

#endif
