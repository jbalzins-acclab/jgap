#ifndef SPARSIFIER_HPP
#define SPARSIFIER_HPP

#include <vector>
#include <Eigen/Dense>

#include "data/BasicDataTypes.hpp"

using namespace std;

namespace jgap {

    class Sparsifier {
    public:
        virtual ~Sparsifier() = default;
        virtual vector<vector<double>> selectSparsePoints(const vector<vector<double>> &allPoints) = 0;
    };
}

#endif //SPARSIFIER_HPP
