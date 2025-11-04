#ifndef UTILS_HPP
#define UTILS_HPP

#include "data/Vector3.hpp"

#include <string>
#include <vector>
#include <Eigen/Dense>
#include <nlohmann/json.hpp>

#include "data/AtomicStructure.hpp"

#define GET_OR_DEFAULT(from, key, defaultValue) !from.contains(key) ? defaultValue : from[key]

using namespace std;

namespace jgap {
    map<string, string> parseHeaderLine(const string &line);
    bool getLine(ifstream &file, string &line);

    vector<AtomicStructure> readXyz(const string& fileName);
    vector<AtomicStructure> readXyz(const string& fileName, double cutoff); // mainly for testing
    void writeXyz(const string& fileName, const vector<AtomicStructure> &structures);

    vector<string> split(const string& s, char delimiter);

    string matrixToString(const Eigen::MatrixXd& mat);
    string vectorToString(const Eigen::VectorXd& vec);
    string vectorToString(const vector<double>& vec);
    string vectorToString(const vector<size_t>& vec);

    nlohmann::json& require(nlohmann::json& j, const string& key);
    const nlohmann::json& require(const nlohmann::json& j, const string& key);

    nlohmann::json requireArray(nlohmann::json &j);
    const nlohmann::json& requireArray(const nlohmann::json& j);

    template<typename Map, typename Key, typename Value>
    auto getOrDefault(const Map& m, const Key& k, const Value& defaultValue) -> decltype(m.at(k)) {
        auto it = m.find(k);
        return it != m.end() ? it->second : defaultValue;
    }

    template <typename Iterator>
    string iteratorToString(Iterator begin, Iterator end) {
        ostringstream oss;
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
