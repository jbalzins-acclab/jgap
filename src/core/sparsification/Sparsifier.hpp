#ifndef SPARSIFIER_HPP
#define SPARSIFIER_HPP

#include <vector>

#include "core/atomic/Descriptor.hpp"

namespace jgap {

    // concept instead?

    template<size_t Dim>
    class Sparsifier {
    public:
        using Ptr = std::shared_ptr<Sparsifier>;

        virtual ~Sparsifier() = default;
        virtual std::vector<Descriptor<Dim>> selectSparsePoints(const std::vector<Descriptor<Dim>> &descriptors) const
            = 0;
    };
}

#endif
