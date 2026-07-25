#ifndef JGAP_UTILS_HPP
#define JGAP_UTILS_HPP

#include <chrono>
#include <cmath>
#include <iosfwd>
#include <map>
#include <ranges>
#include <string>
#include <string_view>
#include <vector>

#include "jgap/core/Real.hpp"

namespace jgap {
    struct MainXYZPropertyNames;
    class Atoms;

#if defined(__GNUC__) || defined(__clang__)
    inline void sincos(float angle, float* sin_val, float* cos_val) { __builtin_sincosf(angle, sin_val, cos_val); }
    inline void sincos(double angle, double* sin_val, double* cos_val) { __builtin_sincos(angle, sin_val, cos_val); }
    inline void sincos(long double angle, long double* sin_val, long double* cos_val) {
        __builtin_sincosl(angle, sin_val, cos_val);
    }
#else
    inline void sincos(float angle, float* sin_val, float* cos_val) {
        *sin_val = std::sin(angle);
        *cos_val = std::cos(angle);
    }
    inline void sincos(double angle, double* sin_val, double* cos_val) {
        *sin_val = std::sin(angle);
        *cos_val = std::cos(angle);
    }
    inline void sincos(long double angle, long double* sin_val, long double* cos_val) {
        *sin_val = std::sin(angle);
        *cos_val = std::cos(angle);
    }
#endif

    bool getLine(std::istream& file, std::string& line);

    std::map<std::string, std::string> parseHeaderLine(const std::string& line);

    std::string uniqueStamp();

    inline double elapsedMillisSince(const std::chrono::steady_clock::time_point& start) {
        auto end = std::chrono::steady_clock::now();
        return std::chrono::duration<double, std::milli>(end - start).count();
    }

    inline std::string formatDuration(double ms) {
        if (ms < 1000.0) {
            return std::to_string(static_cast<int>(ms)) + " ms";
        }
        double sec = ms / 1000.0;
        if (sec < 60.0) {
            return std::to_string(sec) + " s";
        }
        int min = static_cast<int>(sec / 60.0);
        int rem_sec = static_cast<int>(sec) % 60;
        return std::to_string(min) + "m " + std::to_string(rem_sec) + "s";
    }

    constexpr double factorialD(size_t n) {
        double result = 1.0;
        for (size_t i = 2; i <= n; i++) {
            result *= static_cast<Real>(i);
        }
        return result;
    }

    constexpr size_t factorial(size_t n) {
        size_t result = 1;
        for (size_t i = 2; i <= n; i++) {
            result *= i;
        }
        return result;
    }

    double rms(const std::vector<double>&);
    double deviation(const std::vector<double>&);

    std::vector<std::string> split(const std::string& s, char delimiter);
    std::string join(const std::vector<std::string>& s, char delimiter);
    std::string withoutExtension(const std::string& s);

    std::string vectorToString(const std::vector<double>&);
    std::string vectorToString(const std::vector<size_t>&);
    std::string vectorToString(const std::vector<std::string>&);

    template<typename T>
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
        return std::string(name.substr(
            value_start, value_end == std::string_view::npos ? std::string_view::npos : value_end - value_start
        ));
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

    template<typename Iterator>
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

    template<typename InputRange, typename Func>
    auto mapVector(InputRange&& range, Func func) {
        return range | std::views::transform(func) | std::ranges::to<std::vector>();
    }
}

#endif
