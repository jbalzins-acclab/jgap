#ifndef SPARSIFIER_HPP
#define SPARSIFIER_HPP

#include <vector>
#include "jgap/core/atomic/descriptor/Descriptor.hpp"

namespace jgap {

    template<size_t Dim>
    class Sparsifier {
    public:
        virtual ~Sparsifier() = default;
        virtual std::vector<Descriptor<Dim>> selectSparsePoints(
            const std::vector<Descriptor<Dim>>& descriptors) const = 0;
    };
}

#endif
