#ifndef JGAP_MATRIX_HPP
#define JGAP_MATRIX_HPP

#include <cassert>
#include <mdspan>

#include "SharedArray.hpp"

namespace jgap {

    enum class MatrixLayout { RowMajor, ColumnMajor };
    using MatrixLayout::ColumnMajor;
    using MatrixLayout::RowMajor;

    template<MatrixLayout Layout, typename T = double>
    class Matrix {
    public:
        Matrix(const size_t rows, const size_t columns) : rows(rows), columns(columns), memory_space(rows * columns) {}

        Matrix(SharedArray<T> memory_space, const size_t columns) :
            rows(columns > 0 ? memory_space.size() / columns : 0),
            columns(columns),
            memory_space(std::move(memory_space)) {
            assert(columns == 0 || this->memory_space.size() >= columns);
        }

        Matrix(SharedArray<T> memory_space, const size_t rows, const size_t columns) :
            rows(rows), columns(columns), memory_space(std::move(memory_space)) {
            assert(this->memory_space.size() >= rows * columns);
        }

        ~Matrix() = default;

        T& operator()(const size_t i, const size_t j) {
            assert(i < rows && j < columns);
            if constexpr (Layout == MatrixLayout::RowMajor) {
                return memory_space[i * columns + j];
            }
            return memory_space[j * rows + i];
        }

        const T& operator()(const size_t i, const size_t j) const {
            assert(i < rows && j < columns);
            if constexpr (Layout == MatrixLayout::RowMajor) {
                return memory_space[i * columns + j];
            }
            return memory_space[j * rows + i];
        }

        auto mdspan() {
            if constexpr (Layout == MatrixLayout::RowMajor) {
                return std::mdspan<T, std::dextents<size_t, 2>, std::layout_right>(data(), rows, columns);
            } else {
                return std::mdspan<T, std::dextents<size_t, 2>, std::layout_left>(data(), rows, columns);
            }
        }

        auto mdspan() const {
            if constexpr (Layout == MatrixLayout::RowMajor) {
                return std::mdspan<const T, std::dextents<size_t, 2>, std::layout_right>(data(), rows, columns);
            } else {
                return std::mdspan<const T, std::dextents<size_t, 2>, std::layout_left>(data(), rows, columns);
            }
        }

        SharedArray<T>& flatData() { return memory_space; }
        const SharedArray<T>& flatData() const { return memory_space; }

        SharedArray<T>& memorySpace() { return memory_space; }
        const SharedArray<T>& memorySpace() const { return memory_space; }

        T* data() { return memory_space.data(); }
        const T* data() const { return memory_space.data(); }

        size_t nRows() const { return rows; }
        size_t nColumns() const { return columns; }

    private:
        size_t rows{0};
        size_t columns{0};
        SharedArray<T> memory_space;
    };


}

#endif
