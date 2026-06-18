#ifndef JGAP_UTILS_HPP
#define JGAP_UTILS_HPP

#include <string>
#include <vector>
#include <ranges>
#include <cmath>
#include <map>
#include <Eigen/Dense>
#include <string_view>

#include "core/Real.hpp"

namespace jgap {
    struct AtomsPropertyNames;
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

    bool getLine(std::ifstream &file, std::string &line);

    std::vector<Atoms> readAtoms(const std::string& filename);
    std::vector<Atoms> readAtoms(const std::string& filename, const AtomsPropertyNames& names);

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
    constexpr std::string_view typeName() {
#ifdef __clang__
        std::string_view name = __PRETTY_FUNCTION__;
        constexpr std::string_view prefix = "std::string type_name() [T = ";
        constexpr std::string_view suffix = "]";
#elif defined(__GNUC__)
        std::string_view name = __PRETTY_FUNCTION__;
        constexpr std::string_view prefix = "std::string type_name() [with T = ";
        constexpr std::string_view suffix = "]";
#elif defined(_MSC_VER)
        std::string_view name = __FUNCSIG__;
        constexpr std::string_view prefix = "class std::basic_string __cdecl type_name<";
        constexpr std::string_view suffix = ">(void)";
#endif

        // Slice the view safely
        name.remove_prefix(prefix.size());
        name.remove_suffix(suffix.size());

        // Construct and return a real std::string
        return std::string(name);
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