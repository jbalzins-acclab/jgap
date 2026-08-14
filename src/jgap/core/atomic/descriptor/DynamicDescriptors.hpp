#ifndef JGAP_DYNAMICDESCRIPTORS_HPP
#define JGAP_DYNAMICDESCRIPTORS_HPP

#include <cassert>
#include <cstddef>
#include <cstring>
#include <vector>
#include "jgap/core/Real.hpp"
#include "jgap/core/atomic/descriptor/Descriptor.hpp"

namespace jgap {

    /// @brief Store an array of descriptor without knowing their dimensions at compile-time.
    class DynamicDescriptors {
    public:
        DynamicDescriptors(size_t n_descriptors, size_t dim)
            : data(n_descriptors * dim, 0.0),
              n_descriptors(n_descriptors),
              dim(dim) {}

        DynamicDescriptors(std::vector<Real>&& data, size_t dim)
            : data(std::move(data)),
              n_descriptors(dim > 0 ? this->data.size() / dim : 0),
              dim(dim) {
            assert(dim > 0 ? this->data.size() % dim == 0 : this->data.empty());
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
            std::vector<Descriptor<Dim>> result(n_descriptors);
            if (n_descriptors > 0) {
                std::memcpy(result.data(), data.data(), data.size() * sizeof(Real));
            }
            return result;
        }

        template<size_t Dim>
        std::vector<Descriptor<Dim>> releaseAsVector() && {
            assert(Dim == dim);
            std::vector<Descriptor<Dim>> result(n_descriptors);
            if (n_descriptors > 0) {
                std::memcpy(result.data(), data.data(), data.size() * sizeof(Real));
            }
            data.clear();
            n_descriptors = 0;
            return result;
        }

        template<size_t Dim>
        std::vector<Descriptor<Dim>> releaseAsVector() & {
            return toVector<Dim>();
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

        const std::vector<Real>& getData() const {
            return data;
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
