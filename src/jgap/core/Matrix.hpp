#ifndef JGAP_MATRIX_HPP
#define JGAP_MATRIX_HPP

#include <cassert>

#include "SharedArray.hpp"

namespace jgap {

    enum class MatrixLayout { RowMajor, ColumnMajor };
    using MatrixLayout::ColumnMajor;
    using MatrixLayout::RowMajor;

    template<MatrixLayout Layout, typename T = double>
    class Matrix {
    public:
        Matrix() = default;

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

        SharedArray<T>& flatData() { return memory_space; }
        const SharedArray<T>& flatData() const { return memory_space; }

        SharedArray<T>& memorySpace() { return memory_space; }
        const SharedArray<T>& memorySpace() const { return memory_space; }

        T* data() { return memory_space.data(); }
        const T* data() const { return memory_space.data(); }

        size_t nRows() const { return rows; }
        size_t nColumns() const { return columns; }

        void appendRow(const std::vector<T>& row) {
            static_assert(Layout == MatrixLayout::RowMajor, "appendRow currently only supported for RowMajor layout");
            assert(columns == 0 || row.size() == columns);
            if (columns == 0) {
                columns = row.size();
            }
            SharedArray<T> new_space((rows + 1) * columns);
            for (size_t i = 0; i < rows * columns; ++i) {
                new_space[i] = memory_space[i];
            }
            for (size_t j = 0; j < columns; ++j) {
                new_space[rows * columns + j] = row[j];
            }
            memory_space = std::move(new_space);
            ++rows;
        }

        void removeRow(size_t index) {
            static_assert(Layout == MatrixLayout::RowMajor, "removeRow currently only supported for RowMajor layout");
            assert(index < rows);
            SharedArray<T> new_space((rows - 1) * columns);
            size_t dest = 0;
            for (size_t r = 0; r < rows; ++r) {
                if (r == index) continue;
                for (size_t c = 0; c < columns; ++c) {
                    new_space[dest++] = memory_space[r * columns + c];
                }
            }
            memory_space = std::move(new_space);
            --rows;
        }

    private:
        size_t rows{0};
        size_t columns{0};
        SharedArray<T> memory_space;
    };

}

#endif
