#ifndef JGAP_MATRIX_HPP
#define JGAP_MATRIX_HPP

#include <vector>

#include "Real.hpp"
#include "io/log/CurrentLogger.hpp"

namespace jgap {

    /// Simple wrapper for a matrix stored in a single memory block, in a row-major order.
    /// Intended to be used outside the code blocks directly doing the linear algebra
    /// (essentially for \ref GapComponent to provide Sparse-Sparse K_MM matrix blocks).
    class Matrix {
    public:
        Matrix(const size_t rows, const size_t columns)
            : rows(rows), columns(columns) {
            try {
                data = std::vector(rows * columns, Real{});
            } catch (const std::bad_alloc&) {
                JGAP_LOG_AND_THROW("Matrix memory allocation failed");
            }
        }
        ~Matrix() = default;

        Real& operator()(const size_t i, const size_t j) {
            return data[i * columns + j];
        }

        /// \note Emphasizing in name the backend matrix storage layout chosen.
        std::vector<Real>& rowMajorFlatData() { return data; }

        size_t nRows() const { return rows; }
        size_t nColumns() const { return columns; }

    private:
        size_t rows;
        size_t columns;
        std::vector<Real> data;
    };
}

#endif
