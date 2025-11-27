#ifndef UTILS_HPP
#define UTILS_HPP

#include "data/AtomicStructure.hpp"

#include <string>
#include <vector>
#include <Eigen/Dense>
#include <nlohmann/json.hpp>

#define GET_OR_DEFAULT(from, key, defaultValue) !from.contains(key) ? defaultValue : from[key]
#define IO_NOT_INTENDED(func) throw std::logic_error("IO not intended for " + std::string(#func))

namespace jgap {
    std::map<std::string, std::string> parseHeaderLine(const std::string &line);
    bool getLine(std::ifstream &file, std::string &line);

    static std::string uniqueStamp();

    std::vector<AtomicStructure> readXyz(const std::string& fileName);
    std::vector<AtomicStructure> readXyz(const std::string& fileName, double cutoff);
    void writeXyz(const std::string& fileName, const std::vector<AtomicStructure> &structures);
    void writeXyz(std::ofstream& outputStream, const std::vector<AtomicStructure> &structures);

    double rms(const std::vector<double>&);
    double deviation(const std::vector<double>&);

    std::vector<std::string> split(const std::string& s, char delimiter);

    std::string matrixToString(const Eigen::MatrixXd& mat);
    std::string vectorToString(const Eigen::VectorXd& vec);
    std::string vectorToString(const std::vector<double>& vec);
    std::string vectorToString(const std::vector<size_t>& vec);
    std::string vectorToString(const std::vector<std::string>& vec);

    nlohmann::json& requireFull(nlohmann::json& j, const std::string& key, const char* file, int line, const char* function);
    const nlohmann::json& requireFull(const nlohmann::json& j, const std::string& key, const char* file, int line,
                                      const char* function);

    nlohmann::json requireArrayFull(nlohmann::json &j, const char* file, int line, const char* function);
    const nlohmann::json& requireArrayFull(const nlohmann::json& j, const char* file, int line, const char* function);

    #define require(j, key) requireFull((j), (key), __FILE__, __LINE__, __func__)
    #define requireArray(j) requireArrayFull((j), __FILE__, __LINE__, __func__)
    #define optionallySet(val, j, key) if (j.contains(key)) { val = j[key]; }

    template<typename Map, typename Key, typename Value>
    auto getOrDefault(const Map& m, const Key& k, const Value& defaultValue) -> decltype(m.at(k)) {
        auto it = m.find(k);
        return it != m.end() ? it->second : defaultValue;
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
}

#endif
