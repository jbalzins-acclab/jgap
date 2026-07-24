#ifndef JGAP_DESCRIPTORSVIEW_HPP
#define JGAP_DESCRIPTORSVIEW_HPP

#include <cassert>
#include <cstddef>
#include <vector>
#include "jgap/core/Real.hpp"
#include "jgap/core/atomic/descriptor/Descriptor.hpp"

namespace jgap {

    class DescriptorsView {
    public:
        template<size_t Dim>
        DescriptorsView(const std::vector<Descriptor<Dim>>& descriptors)
            : data(descriptors.empty() ? nullptr : descriptors.front().data()),
              n_descriptors(descriptors.size()),
              dim(Dim) {}

        DescriptorsView(const Real* data, size_t n_descriptors, size_t dim)
            : data(data), n_descriptors(n_descriptors), dim(dim) {}

        Real operator()(size_t desc_idx, size_t dim_idx) const {
            assert(desc_idx < n_descriptors);
            assert(dim_idx < dim);
            return data[desc_idx * dim + dim_idx];
        }

        const Real* rawData() const { return data; }
        size_t rawDataSize() const { return n_descriptors * dim; }

        size_t getNumDescriptors() const { return n_descriptors; }
        size_t getDim() const { return dim; }
        bool empty() const { return n_descriptors == 0; }

    private:
        const Real* data;
        size_t n_descriptors;
        size_t dim;
    };

}

#endif
