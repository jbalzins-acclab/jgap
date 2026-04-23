#ifndef JGAP_UTILS_HPP
#define JGAP_UTILS_HPP

#include "data/atomic/AtomicStructure.hpp"

#include <string>
#include <vector>
#include <Eigen/Dense>
#include "../core/DataNode.hpp"

#define GET_OR_DEFAULT(from, key, defaultValue) !from.contains(key) ? defaultValue : from[key]

#define REQUIRE(node, key) requireFull((node), (key), __FILE__, __LINE__, __func__)
#define REQUIRE_ARRAY(node) requireArrayFull((node), __FILE__, __LINE__, __func__)
#define OPTIONALLY_SET(val, node, key) if (node.contains(key)) { val = node[key]; }

namespace jgap {
    class Box;

    std::map<std::string, std::string> parseHeaderLine(const std::string &line);
    bool getLine(std::ifstream &file, std::string &line);

    std::string uniqueStamp();

    double factorial(size_t n);

    std::vector<Box> readXyz(const std::string& file_name);
    void writeXyz(const std::string& fileName, const std::vector<Box> &structures);
    void writeXyz(std::ofstream& outputStream, const std::vector<Box> &structures);

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

    DataNode& requireFull(DataNode& n, const std::string& key, const char* file, int line, const char* function);
    const DataNode& requireFull(const DataNode& n, const std::string& key, const char* file,
                                int line, const char* function);

    DataNode requireArrayFull(DataNode &n, const char* file, int line, const char* function);
    const DataNode& requireArrayFull(const DataNode& n, const char* file, int line, const char* function);

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
}

#endif
