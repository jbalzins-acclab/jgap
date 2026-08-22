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
#include "jgap/core/Vector3.hpp"
#include "jgap/core/atomic/energy/Virials.hpp"
#include "jgap/core/atomic/geometry/Separation.hpp"

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

    double elapsedMillisSince(const std::chrono::steady_clock::time_point& start);

    std::string formatDuration(double ms);

    constexpr Real factorial(size_t n);

    double rms(const std::vector<double>&);
    double deviation(const std::vector<double>&);

    std::vector<std::string> split(const std::string& s, char delimiter);
    std::string join(const std::vector<std::string>& s, char delimiter);
    std::string withoutExtension(const std::string& s);

    std::string vectorToString(const std::vector<double>&);
    std::string vectorToString(const std::vector<size_t>&);
    std::string vectorToString(const std::vector<std::string>&);

    inline void accumulatePairDerivatives(Vector3& f0, Vector3& f1, Virials& virials,
                                          const Real dE_dr, const Separation& r01) {
        virials += r01.virials() * dE_dr;

        Vector3 f10 = r01.direction * dE_dr;
        f0 += f10;
        f1 -= f10;
    }

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
