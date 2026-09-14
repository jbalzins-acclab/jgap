#ifndef JGAP_BLOCKINCREMENTALQRACCUMULATOR_HPP
#define JGAP_BLOCKINCREMENTALQRACCUMULATOR_HPP

#include <mutex>
#include <optional>
#include <vector>

#include "jgap/core/Matrix.hpp"
#include "jgap/core/SharedArray.hpp"


namespace jgap::linalg {

    class BlockIncrementalQRAccumulator {
    public:
        static SharedArray<double> allocateMatrixMemory(size_t n_cols, double approx_ram_limit_gb);

        BlockIncrementalQRAccumulator(
            SharedArray<double> A,
            SharedArray<double> b,
            size_t n_cols,
            size_t n_rows_filled
        );

        BlockIncrementalQRAccumulator(
            const Matrix<ColumnMajor>& A,
            SharedArray<double> b,
            size_t n_rows_filled
        );

        size_t nIncrementBlockRows() const { return n_increment_block_rows; }

        void appendBlock(const Matrix<ColumnMajor>& A_block, const std::vector<double>& b_chunk);
        std::vector<double> solve();

        void finalize();

        SharedArray<double> getMatrixMemorySpace() const { return A.flatData(); }

    private:
        size_t n_rows_filled;
        size_t n_cols;
        size_t n_increment_block_rows;
        Matrix<ColumnMajor> A;
        SharedArray<double> b;
        std::mutex mutex;

        static size_t calculateMaxIncrementBlockRows(size_t n_cols, double approx_ram_limit_gb);
        void flushFullBlock();
    };
}

#endif
