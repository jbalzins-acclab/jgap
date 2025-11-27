#ifndef SPARSIFIER_HPP
#define SPARSIFIER_HPP

#include <vector>
#include <Eigen/Dense>

#include "data/Vector3.hpp"

namespace jgap {

    using KernelParams = nlohmann::json;

    class Sparsifier {
    public:
        virtual ~Sparsifier() = default;
        virtual std::vector<KernelParams> selectSparsePoints(const std::vector<std::vector<double>> &allPoints) = 0;
    };
}

#endif
