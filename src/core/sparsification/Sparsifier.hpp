#ifndef SPARSIFIER_HPP
#define SPARSIFIER_HPP

#include <vector>

#include "data/DataNode.hpp"

namespace jgap {

    template<size_t N_DIMENSIONS, size_t N_ATOMS>
    class Sparsifier {
    public:
        virtual ~Sparsifier() = default;
        virtual std::vector<Descriptor<N_DIMENSIONS, N_ATOMS>> selectSparsePoints(
            const std::vector<Descriptor<N_DIMENSIONS, N_ATOMS>> &descriptors) = 0;
    };
}

#endif
