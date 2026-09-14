#ifndef JGAP_SHAREDARRAY_HPP
#define JGAP_SHAREDARRAY_HPP

#include <cassert>
#include <memory>

#include "jgap/core/io/log/CurrentLogger.hpp"

namespace jgap {

    template<typename T>
    struct SharedArray {
        std::shared_ptr<T[]> ptr{nullptr};
        size_t n_elements{0};

        SharedArray() = default;

        explicit SharedArray(const size_t n_elements) : n_elements(n_elements) {
            try {
                ptr = std::make_shared<T[]>(n_elements);
            } catch (const std::bad_alloc&) {
                JGAP_LOG_AND_THROW("SharedArray memory allocation failed");
            }
        }

        SharedArray(std::shared_ptr<T[]> ptr, const size_t n_elements) :
            ptr(std::move(ptr)), n_elements(n_elements) {}

        SharedArray subspace(const size_t starting_element) const {
            assert(starting_element <= n_elements);
            const size_t sub_count = n_elements - starting_element;
            std::shared_ptr<T[]> sub_ptr(ptr, ptr.get() + starting_element);
            return SharedArray(std::move(sub_ptr), sub_count);
        }

        SharedArray subspace(const size_t starting_element, const size_t last_element) const {
            assert(starting_element <= last_element);
            assert(last_element <= n_elements);
            const size_t sub_count = last_element - starting_element;
            std::shared_ptr<T[]> sub_ptr(ptr, ptr.get() + starting_element);
            return SharedArray(std::move(sub_ptr), sub_count);
        }

        T* get() const { return ptr.get(); }
        T* data() const { return ptr.get(); }
        std::shared_ptr<T[]> sharedPtr() const { return ptr; }

        T& operator[](const size_t i) {
            assert(i < n_elements);
            return ptr[i];
        }

        const T& operator[](const size_t i) const {
            assert(i < n_elements);
            return ptr[i];
        }

        T* begin() const { return ptr.get(); }
        T* end() const { return ptr.get() + n_elements; }

        size_t size() const { return n_elements; }

        void fill(const T& value) {
            for (size_t i = 0; i < n_elements; i++) {
                ptr[i] = value;
            }
        }
    };

}

#endif
