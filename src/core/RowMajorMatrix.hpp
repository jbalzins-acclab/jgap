#ifndef JGAP_MATRIXBLOCK_HPP
#define JGAP_MATRIXBLOCK_HPP

#include <vector>

#include "Real.hpp"
#include "io/log/CurrentLogger.hpp"

namespace jgap {

    /**
     * Simple wrapper for a matrix stored in a single memory block, in a row-major order.
     * Intended to be used outside the code blocks directly doing the linear algebra
     * (essentially for @link GapComponent to provide Sparse-Sparse K_MM matrix blocks).
     */
    class RowMajorMatrix {
    public:
        RowMajorMatrix(const size_t rows, const size_t columns)
            : rows(rows), columns(columns) {
            try {
                data = std::vector(rows * columns, Real{});
            } catch (const std::bad_alloc&) {
                JGAP_LOG_AND_THROW("Matrix memory allocation failed");
            }
        }
        ~RowMajorMatrix() = default;

        Real& operator()(const size_t i, const size_t j) {
            return data[i * columns + j];
        }

        std::vector<Real>& rawData() { return data; }

        size_t nRows() const { return rows; }
        size_t nColumns() const { return columns; }
        std::pair<size_t, size_t> blockSize() const { return{rows, columns}; }

    private:
        size_t rows;
        size_t columns;
        std::vector<Real> data;
    };
}

#endif
