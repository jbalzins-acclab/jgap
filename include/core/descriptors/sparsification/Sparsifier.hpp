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
        virtual vector<Eigen::VectorXd> selectSparsePoints(const vector<Eigen::VectorXd> &allPoints) = 0;
    };
}

#endif //SPARSIFIER_HPP
