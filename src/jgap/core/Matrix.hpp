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

        void makeContiguous() {
            if constexpr (Layout == MatrixLayout::ColumnMajor) {
                const size_t N = columns;
                assert(rows >= N);
                for (size_t j = 0; j < N; ++j) {
                    for (size_t i = 0; i < N; ++i) {
                        T val = (i > j) ? T{0} : memory_space[j * rows + i];
                        memory_space[j * N + i] = val;
                    }
                }
            } else {
                const size_t N = rows;
                assert(columns >= N);
                for (size_t i = 0; i < N; ++i) {
                    for (size_t j = 0; j < N; ++j) {
                        T val = (i > j) ? T{0} : memory_space[i * columns + j];
                        memory_space[i * N + j] = val;
                    }
                }
            }
        }

        void makeContagious() { makeContiguous(); }

        static void unstack(
            const Matrix<Layout, T>& source,
            Matrix<Layout, T>& target,
            const size_t target_start_row,
            const size_t target_start_col
        ) {
            if constexpr (Layout == MatrixLayout::ColumnMajor) {
                const size_t N = source.nColumns();
                assert(source.nRows() == N);
                assert(target_start_row + N <= target.nRows());
                assert(target_start_col + N <= target.nColumns());

                for (size_t j_idx = N; j_idx > 0; --j_idx) {
                    const size_t j = j_idx - 1;
                    const size_t target_c = target_start_col + j;

                    for (size_t r = target.nRows(); r > target_start_row + N; --r) {
                        target(r - 1, target_c) = T{0};
                    }

                    for (size_t i_idx = N; i_idx > 0; --i_idx) {
                        const size_t i = i_idx - 1;
                        const T val = (i > j) ? T{0} : source(i, j);
                        target(target_start_row + i, target_c) = val;
                    }

                    for (size_t r = target_start_row; r > 0; --r) {
                        target(r - 1, target_c) = T{0};
                    }
                }
            } else {
                const size_t N = source.nRows();
                assert(source.nColumns() == N);
                assert(target_start_row + N <= target.nRows());
                assert(target_start_col + N <= target.nColumns());

                for (size_t i_idx = N; i_idx > 0; --i_idx) {
                    const size_t i = i_idx - 1;
                    const size_t target_r = target_start_row + i;

                    for (size_t c = target.nColumns(); c > target_start_col + N; --c) {
                        target(target_r, c - 1) = T{0};
                    }

                    for (size_t j_idx = N; j_idx > 0; --j_idx) {
                        const size_t j = j_idx - 1;
                        const T val = (i > j) ? T{0} : source(i, j);
                        target(target_r, target_start_col + j) = val;
                    }

                    for (size_t c = target_start_col; c > 0; --c) {
                        target(target_r, c - 1) = T{0};
                    }
                }
            }
        }

    private:
        size_t rows{0};
        size_t columns{0};
        SharedArray<T> memory_space;
    };


}

#endif
