#ifndef JGAP_DYNAMICDESCRIPTORS_HPP
#define JGAP_DYNAMICDESCRIPTORS_HPP

#include <vector>
#include <cstddef>
#include <cassert>
#include "core/Real.hpp"
#include "core/atomic/descriptor/Descriptor.hpp"

namespace jgap {

    class DynamicDescriptors {
    public:
        DynamicDescriptors(size_t n_descriptors, size_t dim)
            : data(n_descriptors * dim, 0.0),
              n_descriptors(n_descriptors),
              dim(dim) {}

        DynamicDescriptors(std::vector<Real>&& data, size_t dim)
            : data(std::move(data)),
              n_descriptors(this->data.size() / dim),
              dim(dim) {
            assert(this->data.size() % dim == 0);
        }

        Real& operator()(size_t desc_idx, size_t dim_idx) {
            assert(desc_idx < n_descriptors);
            assert(dim_idx < dim);
            return data[desc_idx * dim + dim_idx];
        }

        Real operator()(size_t desc_idx, size_t dim_idx) const {
            assert(desc_idx < n_descriptors);
            assert(dim_idx < dim);
            return data[desc_idx * dim + dim_idx];
        }

        template<size_t Dim>
        std::vector<Descriptor<Dim>> toVector() const {
            assert(Dim == dim);
            std::vector<Descriptor<Dim>> result;
            result.reserve(n_descriptors);
            for (size_t i = 0; i < n_descriptors; ++i) {
                Descriptor<Dim> d;
                for (size_t j = 0; j < Dim; ++j) {
                    d.value[j] = (*this)(i, j);
                }
                result.push_back(d);
            }
            return result;
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

        void pushBack(const std::vector<Real>& point) {
            assert(point.size() == dim);
            data.insert(data.end(), point.begin(), point.end());
            n_descriptors++;
        }

        void pushBackRaw(const Real* point) {
            data.insert(data.end(), point, point + dim);
            n_descriptors++;
        }

    private:
        std::vector<Real> data;
        size_t n_descriptors;
        size_t dim;
    };

}

#endif
