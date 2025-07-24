#ifndef UTILS_HPP
#define UTILS_HPP

#include "data/BasicDataTypes.hpp"

#include <string>
#include <vector>
#include <Eigen/Dense>

using namespace std;

namespace jgap {
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
}

#endif
