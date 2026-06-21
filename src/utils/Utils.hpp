#ifndef JGAP_UTILS_HPP
#define JGAP_UTILS_HPP

#include <string>
#include <vector>
#include <ranges>
#include <cmath>
#include <chrono>
#include <map>
#include <iosfwd>
#include <Eigen/Dense>
#include <string_view>

#include "core/Real.hpp"

namespace jgap {
    struct MainXYZPropertyNames;
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

    bool getLine(std::istream &file, std::string &line);

    std::map<std::string, std::string> parseHeaderLine(const std::string &line);

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

    template <typename T>
    std::string typeName() {
#if defined(__clang__) || defined(__GNUC__)
        std::string_view name = __PRETTY_FUNCTION__;
        // Extract the substring between "T = " and the closing ']'/';'.
        const auto start = name.find("T = ");
        if (start == std::string_view::npos) {
            return std::string(name);
        }
        const auto value_start = start + 4;
        const auto value_end = name.find_first_of(";]", value_start);
        return std::string(name.substr(value_start,
            value_end == std::string_view::npos ? std::string_view::npos : value_end - value_start));
#elif defined(_MSC_VER)
        std::string_view name = __FUNCSIG__;
        constexpr std::string_view prefix = "typeName<";
        constexpr std::string_view suffix = ">(void)";
        const auto start = name.find(prefix);
        if (start == std::string_view::npos) {
            return std::string(name);
        }
        const auto value_start = start + prefix.size();
        const auto value_end = name.rfind(suffix);
        return std::string(name.substr(value_start, value_end - value_start));
#else
        return std::string(__func__);
#endif
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