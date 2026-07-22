#ifndef JGAP_DESCRIPTORSVIEW_HPP
#define JGAP_DESCRIPTORSVIEW_HPP

#include <vector>
#include <cstddef>
#include <cassert>
#include "core/Real.hpp"
#include "core/atomic/descriptor/Descriptor.hpp"

namespace jgap {

    class DescriptorsView {
    public:
        template<size_t Dim>
        explicit DescriptorsView(const std::vector<Descriptor<Dim>>& descriptors)
            : data(descriptors.empty() ? nullptr : descriptors.front().data()),
              n_descriptors(descriptors.size()),
              dim(Dim) {}

        Real operator()(size_t desc_idx, size_t dim_idx) const {
            assert(desc_idx < n_descriptors);
            assert(dim_idx < dim);
            return data[desc_idx * dim + dim_idx];
        }

        size_t getNumDescriptors() const {
            return n_descriptors;
        }

        size_t getDim() const {
            return dim;
        }

        bool empty() const {
            return n_descriptors == 0;
        }

    private:
        const Real* data;
        size_t n_descriptors;
        size_t dim;
    };

}

#endif
