#ifndef UTILS_HPP
#define UTILS_HPP

#include "data/BasicDataTypes.hpp"

#include <string>
#include <vector>
#include <Eigen/Dense>

#define GET_OR_DEFAULT(from, key, defaultValue) !from.contains(key) ? defaultValue : from[key]

using namespace std;

namespace jgap {
    map<string, string> parseHeaderLine(const string &line);
    bool getLine(ifstream &file, string &line);
    vector<AtomicStructure> readXyz(const string& fileName);
    vector<AtomicStructure> readXyz(const string& fileName, double cutoff); // mainly for testing
    void writeXyz(const string& fileName, const vector<AtomicStructure> &);
    vector<string> split(const string& s, char delimiter);
    void saveArray(const vector<double>& data, const string& filename);
    vector<double> loadArray(const string& filename);
    string matrixToString(const Eigen::MatrixXd& mat);
    string vectorToString(const Eigen::VectorXd& vec);
    string vectorToString(const vector<double>& vec);
    string vectorToString(const vector<size_t>& vec);

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
