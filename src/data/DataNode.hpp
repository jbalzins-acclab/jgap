#ifndef JGAP_DATANODE_HPP
#define JGAP_DATANODE_HPP

#include <cassert>
#include <variant>
#include <vector>
#include <map>
#include <string>
#include <stdexcept>
#include <initializer_list>
#include <type_traits>
#include <utility>
#include <yaml-cpp/yaml.h>

#include "io/log/CurrentLogger.hpp"

namespace jgap {

    struct DataNode {
        enum class Type { BOOL, INT, DOUBLE, STRING, ARRAY, OBJECT };

        static DataNode object() { return DataNode(std::map<std::string, DataNode>{}); }
        static DataNode array()  { return DataNode(std::vector<DataNode>{}); }

        Type type;
        std::variant<
            bool,
            int64_t,
            double,
            std::string,
            std::vector<DataNode>,
            std::map<std::string, DataNode>
        > value;

        DataNode() { value = {}; type = Type::OBJECT; };
        DataNode(const DataNode& cpy) { value = cpy.value; type = cpy.type; };

        // --- Constructors for convenience ---
        DataNode(bool b) : type(Type::BOOL), value(b) {}
        DataNode(int64_t i) : type(Type::INT), value(i) {}
        DataNode(int i) : type(Type::INT), value(int64_t(i)) {}
        DataNode(double d) : type(Type::DOUBLE), value(d) {}
        DataNode(const char* s) : type(Type::STRING), value(std::string(s)) {}
        DataNode(const std::string& s) : type(Type::STRING), value(s) {}
        DataNode(const std::vector<DataNode>& v) : type(Type::ARRAY), value(v) {}
        DataNode(const std::map<std::string, DataNode>& m) : type(Type::OBJECT), value(m) {}

        // Construct array from initializer_list of scalar-like or std::string values
        template<typename T,
                 typename = std::enable_if_t<
                     std::is_arithmetic_v<T> || std::is_same_v<std::decay_t<T>, std::string>
                 >
        >
        DataNode(std::initializer_list<T> init) : type(Type::ARRAY) {
            std::vector<DataNode> v;
            v.reserve(init.size());
            for (const auto& x : init) {
                if constexpr (std::is_same_v<std::decay_t<T>, std::string>) {
                    v.emplace_back(DataNode{x});
                } else if constexpr (std::is_floating_point_v<T>) {
                    v.emplace_back(DataNode{static_cast<double>(x)});
                } else if constexpr (std::is_integral_v<T>) {
                    v.emplace_back(DataNode{static_cast<int64_t>(x)});
                } else {
                    static_assert(sizeof(T) == 0, "Unsupported initializer_list element type for DataNode array");
                }
            }
            value = std::move(v);
        }

        // Construct object from initializer_list of std::pair-like { key, value }
        DataNode(std::initializer_list<std::pair<std::string, DataNode>> init)
            : type(Type::OBJECT) {
            std::map<std::string, DataNode> m;
            for (const auto& kv : init) {
                m[kv.first] = kv.second; // last wins on duplicate keys
            }
            value = std::move(m);
        }

        template<typename K, typename V, typename = std::enable_if_t<std::is_convertible_v<K, std::string>>>
        DataNode(std::initializer_list<std::pair<K, V>> init)
            : type(Type::OBJECT) {
            std::map<std::string, DataNode> m;
            for (const auto& kv : init) {
                m[static_cast<std::string>(kv.first)] = DataNode(kv.second);
            }
            value = std::move(m);
        }

        // Convenience constructors from common vectors
        DataNode(const std::vector<std::string>& v) : type(Type::ARRAY) {
            std::vector<DataNode> arr; arr.reserve(v.size());
            for (const auto& x : v) arr.emplace_back(DataNode{x});
            value = std::move(arr);
        }
        // Generic std::vector of arithmetic (except bool) maps to array
        template<typename T,
                 typename = std::enable_if_t<std::is_arithmetic_v<T> && !std::is_same_v<T, bool>>>
        DataNode(const std::vector<T>& v) : type(Type::ARRAY) {
            std::vector<DataNode> arr; arr.reserve(v.size());
            for (const auto& x : v) {
                if constexpr (std::is_floating_point_v<T>) arr.emplace_back(DataNode{static_cast<double>(x)});
                else arr.emplace_back(DataNode{static_cast<int64_t>(x)});
            }
            value = std::move(arr);
        }

        // --- Add convenient operator[] overloads that accept string_view / char* ---
        // Non-const (creates element when missing)
        DataNode& operator[](std::string key) {
            if (type != Type::OBJECT)
                throw std::runtime_error("Node is not an object");
            auto& m = std::get<std::map<std::string, DataNode>>(value);
            return m[key]; // copy into key (std::map owns the std::string)
        }

        // Const version (throws if missing)
        const DataNode& operator[](std::string key) const {
            if (type != Type::OBJECT)
                throw std::runtime_error("Node is not an object");
            const auto& m = std::get<std::map<std::string, DataNode>>(value);
            auto it = m.find(std::string(key));
            if (it == m.end()) throw std::runtime_error("Key not found: " + std::string(key));
            return it->second;
        }

        DataNode& operator[](size_t idx) {
            if (type != Type::ARRAY)
                throw std::runtime_error("Node is not an array");
            return std::get<std::vector<DataNode>>(value)[idx];
        }

        bool contains(const std::string& key) const {
            assert (type == Type::OBJECT && "Node is not an object");
            const auto& m = std::get<std::map<std::string, DataNode>>(value);
            return m.find(key) != m.end();
        }
        const DataNode& operator[](size_t idx) const {
            assert (type == Type::ARRAY && "Node is not an array");
            return std::get<std::vector<DataNode>>(value)[idx];
        }

        void pushBack(const DataNode& elem) {
            if (type != Type::ARRAY)
                throw std::runtime_error("pushBack() on non-array DataNode");

            std::get<std::vector<DataNode>>(value).push_back(elem);
        }

        DataNode& back() {
            if (type != Type::ARRAY)
                throw std::runtime_error("back() on non-array DataNode");

            auto& v = std::get<std::vector<DataNode>>(value);
            if (v.empty())
                throw std::runtime_error("back() on empty DataNode array");

            return v.back();
        }

        const DataNode& back() const {
            if (type != Type::ARRAY)
                throw std::runtime_error("back() on non-array DataNode");

            const auto& v = std::get<std::vector<DataNode>>(value);
            if (v.empty())
                throw std::runtime_error("back() on empty DataNode array");

            return v.back();
        }

        // --- Size helpers ---
        size_t size() const {
            if (type == Type::ARRAY) return std::get<std::vector<DataNode>>(value).size();
            if (type == Type::OBJECT) return std::get<std::map<std::string, DataNode>>(value).size();
            return 0;
        }

        template<typename T>
        T as() const {
            return std::get<T>(value);
        }

        auto&       asObject()       { return std::get<std::map<std::string, DataNode>>(value); }
        const auto& asObject() const { return std::get<std::map<std::string, DataNode>>(value); }
        auto&       asArray()        { return std::get<std::vector<DataNode>>(value); }
        const auto& asArray()  const { return std::get<std::vector<DataNode>>(value); }
        std::string asString() const { return std::get<std::string>(value); }
        int64_t     asInt()    const { return std::get<int64_t>(value); }
        double      asDouble() const { return std::get<double>(value); }
        bool        asBool()   const { return std::get<bool>(value); }

        operator std::map<std::string, DataNode>() const { return std::get<std::map<std::string, DataNode>>(value); }
        operator std::vector<DataNode>() const { return std::get<std::vector<DataNode>>(value); }
        operator std::string() const { return std::get<std::string>(value); }
        operator int64_t()     const { return std::get<int64_t>(value); }
        operator double()      const { return std::get<double>(value); }
        operator bool()        const { return std::get<bool>(value); }

        template<typename T>
        static T valueIfMatchesType(const DataNode& node, const T& _typeOf) {
            // First: if the node is null, treat like missing → return default

            try {
                // --- IDENTICAL TYPE DIRECT MATCH ---
                if constexpr (std::is_same_v<T, std::string>) {
                    if (node.type == Type::STRING)
                        return std::get<std::string>(node.value);
                }
                if constexpr (std::is_same_v<T, bool>) {
                    if (node.type == Type::BOOL)
                        return std::get<bool>(node.value);
                }
                if constexpr (std::is_integral_v<T> && !std::is_same_v<T,bool>) {
                    if (node.type == Type::INT)
                        return static_cast<T>(std::get<int64_t>(node.value));
                    if (node.type == Type::DOUBLE)
                        return static_cast<T>(std::get<double>(node.value));
                }
                if constexpr (std::is_floating_point_v<T>) {
                    if (node.type == Type::DOUBLE)
                        return static_cast<T>(std::get<double>(node.value));
                    if (node.type == Type::INT)
                        return static_cast<T>(std::get<int64_t>(node.value));
                }

                // --- SCALAR → STRING ---
                if constexpr (std::is_same_v<T, std::string>) {
                    switch (node.type) {
                        case Type::BOOL:
                            return std::get<bool>(node.value) ? "true" : "false";
                        case Type::INT:
                            return std::to_string(std::get<int64_t>(node.value));
                        case Type::DOUBLE:
                            return std::to_string(std::get<double>(node.value));
                        default: break;
                    }
                }

                // --- STRING → BOOL/INT/DOUBLE ---
                if (node.type == Type::STRING) {
                    const auto& s = std::get<std::string>(node.value);

                    if constexpr (std::is_same_v<T,bool>) {
                        if (s == "true") return true;
                        if (s == "false") return false;
                    }
                    if constexpr (std::is_integral_v<T>) {
                        try { return static_cast<T>(std::stoll(s)); } catch (...) {}
                    }
                    if constexpr (std::is_floating_point_v<T>) {
                        try { return static_cast<T>(std::stod(s)); } catch (...) {}
                    }
                }

            } catch (...) {
                throw std::runtime_error("Failed to convert DataNode value to target type.");
            }

            // If we reach here: conversion impossible → THROW
            throw std::runtime_error(
                "Cannot convert DataNode(" +
                std::string(typeName(node.type)) +   // helper function below
                ") to requested type."
            );
        }

        // Optional helper to turn enum into std::string for debugging:
        static const char* typeName(Type t) {
            switch (t) {
                case Type::BOOL:   return "BOOL";
                case Type::INT:    return "INT";
                case Type::DOUBLE: return "DOUBLE";
                case Type::STRING: return "STRING";
                case Type::ARRAY:  return "ARRAY";
                case Type::OBJECT: return "OBJECT";
            }
            return "UNKNOWN";
        }

        template<typename T>
        T getOrDefault(const std::string& key, const T& defaultValue) const {
            if (type != Type::OBJECT) return defaultValue;
            const auto& m = std::get<std::map<std::string, DataNode>>(value);
            auto it = m.find(key);
            if (it == m.end()) return defaultValue;
            return valueIfMatchesType(it->second, defaultValue);
        }

        // --- Utility methods for range-based loops ---

        // Check if empty
        bool empty() const {
            if (type == Type::ARRAY) return std::get<std::vector<DataNode>>(value).empty();
            if (type == Type::OBJECT) return std::get<std::map<std::string, DataNode>>(value).empty();
            return true;
        }

        // Convert to std::string representation (JSON-like format)
        std::string toString() const {
            switch (type) {
                case Type::BOOL:
                    return std::get<bool>(value) ? "true" : "false";
                case Type::INT:
                    return std::to_string(std::get<int64_t>(value));
                case Type::DOUBLE:
                    return std::to_string(std::get<double>(value));
                case Type::STRING:
                    return "\"" + std::get<std::string>(value) + "\"";
                case Type::ARRAY: {
                    const auto& arr = std::get<std::vector<DataNode>>(value);
                    std::string result = "[";
                    for (size_t i = 0; i < arr.size(); ++i) {
                        if (i > 0) result += ", ";
                        result += arr[i].toString();
                    }
                    result += "]";
                    return result;
                }
                case Type::OBJECT: {
                    const auto& obj = std::get<std::map<std::string, DataNode>>(value);
                    std::string result = "{";
                    bool first = true;
                    for (const auto& [key, val] : obj) {
                        if (!first) result += ", ";
                        first = false;
                        result += key + ": " + val.toString();
                    }
                    result += "}";
                    return result;
                }
            }
            JGAP_LOG_AND_THROW("Failed to serialize as string");
        }
    };
}

#endif