#ifndef JGAP_UTILS_HPP
#define JGAP_UTILS_HPP

#include <string>
#include <vector>
#include <ranges>
#include <cmath>
#include <Eigen/Dense>

#include "core/Real.hpp"

#define GET_OR_DEFAULT(from, key, defaultValue) !from.contains(key) ? defaultValue : from[key]

#define REQUIRE(node, key) requireFull((node), (key), __FILE__, __LINE__, __func__)
#define REQUIRE_ARRAY(node) requireArrayFull((node), __FILE__, __LINE__, __func__)
#define OPTIONALLY_SET(val, node, key) if (node.contains(key)) { val = node[key]; }

namespace jgap {
    class Atoms;

    // Portable sincos implementation
    #if defined(__GNUC__) || defined(__clang__)
        inline void sincos(Real angle, Real* sin_val, Real* cos_val) {
            __builtin_sincos(angle, sin_val, cos_val);
        }
    #else
        inline void sincos(Real angle, Real* sin_val, Real* cos_val) {
            *sin_val = std::sin(angle);
            *cos_val = std::cos(angle);
        }
    #endif

    std::vector<Atoms> readAtoms(const std::string& filename);

    std::string uniqueStamp();

    double factorial(size_t n);

    double rms(const std::vector<double>&);
    double deviation(const std::vector<double>&);

    std::vector<std::string> split(const std::string& s, char delimiter);
    std::string join(const std::vector<std::string>& s, char delimiter);
    std::string withoutExtension(const std::string& s);

    std::string matrixToString(const Eigen::MatrixXd&);
    std::string vectorToString(const Eigen::VectorXd&);
    std::string vectorToString(const std::vector<double>&);
    std::string vectorToString(const std::vector<size_t>&);
    std::string vectorToString(const std::vector<std::string>&);

    /*
    DataNode& requireFull(DataNode& n, const std::string& key, const char* file, int line, const char* function);
    const DataNode& requireFull(const DataNode& n, const std::string& key, const char* file,
                                int line, const char* function);

    DataNode requireArrayFull(DataNode &n, const char* file, int line, const char* function);
    const DataNode& requireArrayFull(const DataNode& n, const char* file, int line, const char* function);
    */
    template<typename Map, typename Key, typename Value>
    auto getOrDefault(const Map& m, const Key& k, const Value& default_value) -> decltype(m.at(k)) {
        auto it = m.find(k);
        return it != m.end() ? it->second : default_value;
    }

    template <typename Iterator>
    std::string iteratorToString(Iterator begin, Iterator end) {
        std::ostringstream oss;
        oss << "[";

        if (begin != end) {
            oss << *begin;
            ++begin;
        }

        while (begin != end) {
            oss << ", " << *begin;
            ++begin;
        }

        oss << "]";
        return oss.str();
    }

    template <typename InputRange, typename Func>
    auto mapVector(InputRange&& range, Func func) {
        return range | std::views::transform(func) | std::ranges::to<std::vector>();
    }
}

#endif