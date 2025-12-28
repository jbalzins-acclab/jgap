#ifndef SPARSIFIER_HPP
#define SPARSIFIER_HPP

#include <vector>

#include "data/DataNode.hpp"

namespace jgap {

    using KernelParams = DataNode;

    class Sparsifier {
    public:
        virtual ~Sparsifier() = default;
        virtual std::vector<KernelParams> selectSparsePoints(const std::vector<std::vector<double>> &allPoints) = 0;
    };
}

#endif
