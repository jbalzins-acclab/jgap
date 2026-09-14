#include "jgap/core/linalg/BlockIncrementalQRAccumulator.hpp"

#include <Eigen/Dense>
#include <Eigen/QR>
#include <cassert>
#include <vector>

#include "jgap/core/io/log/CurrentLogger.hpp"

namespace jgap::linalg {

    using EigenMatrixColMajor = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>;
    using EigenVector = Eigen::Matrix<double, Eigen::Dynamic, 1>;

    size_t BlockIncrementalQRAccumulator::calculateMaxIncrementBlockRows(const size_t n_cols, const double approx_ram_limit_gb) {
        const double max_bytes = approx_ram_limit_gb * 1024.0 * 1024.0 * 1024.0;
        const double bytes_per_row = static_cast<double>(n_cols) * sizeof(double);
        const double mm_bytes = static_cast<double>(n_cols) * bytes_per_row;
        const double min_required_bytes = mm_bytes + bytes_per_row;

        if (max_bytes < min_required_bytes) {
            JGAP_LOG_WARN(
                "Requested RAM limit ({:.2f} GB) is less than the minimum required for full (M+1)xM in-place QR "
                "workspace "
                "({:.2f} GB for M = {}). Using minimal chunk size of 1 row.",
                approx_ram_limit_gb,
                min_required_bytes / (1024.0 * 1024.0 * 1024.0),
                n_cols
            );
            return 1;
        }

        const double remaining_bytes = max_bytes - mm_bytes;
        const auto max_b = static_cast<size_t>(remaining_bytes / bytes_per_row);
        return std::max<size_t>(1, max_b);
    }

    SharedArray<double> BlockIncrementalQRAccumulator::allocateMatrixMemory(size_t n_cols, double approx_ram_limit_gb) {
        const size_t total_rows = n_cols + calculateMaxIncrementBlockRows(n_cols, approx_ram_limit_gb);
        const size_t workspace_bytes = total_rows * n_cols;
        return SharedArray<double>(workspace_bytes);
    }

    BlockIncrementalQRAccumulator::BlockIncrementalQRAccumulator(
        SharedArray<double> A, SharedArray<double> b, const size_t n_cols, size_t n_rows_filled
    ) :
        n_cols(n_cols),
        n_increment_block_rows(n_cols > 0 && A.size() >= n_cols * n_cols ? (A.size() / n_cols) - n_cols : 0),
        n_rows_filled(n_rows_filled),
        A(std::move(A), n_cols),
        b(std::move(b)) {
        assert(n_cols > 0);
        assert(this->A.nRows() >= n_cols);
        assert(b.size() >= this->A.nRows());
    }

    BlockIncrementalQRAccumulator::BlockIncrementalQRAccumulator(
        const Matrix<ColumnMajor>& A, SharedArray<double> b, size_t n_rows_filled
        ) :
        n_cols(A.nColumns()),
        n_increment_block_rows(A.nRows() - A.nColumns()),
        n_rows_filled(n_rows_filled),
        A(A),
        b(std::move(b)) {
        assert(n_cols > 0);
        assert(this->A.nRows() >= n_cols);
        assert(b.size() >= this->A.nRows());
    }

    void BlockIncrementalQRAccumulator::flushFullBlock() {
        const size_t total_rows = n_cols + n_increment_block_rows;
        assert(n_rows_filled == total_rows);

        Eigen::Map<EigenMatrixColMajor> A_eigen(A.data(), total_rows, n_cols);
        Eigen::Map<EigenVector> b_eigen(b.data(), total_rows);

        EigenVector householder_coeffs(n_cols);

        // In-place blocked Householder QR directly inside A_map (0 extra matrix allocations!)
        Eigen::internal::householder_qr_inplace_blocked<Eigen::Map<EigenMatrixColMajor>, EigenVector>::run(
            A_eigen, householder_coeffs
        );

        // Apply Q^T in-place to the target vector using the implicit Householder sequence
        Eigen::HouseholderSequence<Eigen::Map<EigenMatrixColMajor>, EigenVector> householder_Q(
            A_eigen, householder_coeffs
        );
        b_eigen = householder_Q.transpose() * b_eigen;

        // Zero out the strictly lower-triangular part of the top n_cols x n_cols block of R
        A_eigen.topRows(n_cols).template triangularView<Eigen::StrictlyLower>().setZero();

        n_rows_filled = n_cols;
    }

    void BlockIncrementalQRAccumulator::finalize() {
        if (n_rows_filled < n_cols) {
            JGAP_LOG_AND_THROW("Not enough rows provided to perform QR decomposition.");
        }
        if (n_rows_filled == n_cols) {
            return;
        }
        const size_t total_rows = n_cols + n_increment_block_rows;

        // Use explicit OuterStride for the n_rows_filled x n_cols submatrix to guarantee
        // optimal vectorization and avoid complex nested template types in Householder QR.
        using StridedMatrixMap = Eigen::Map<EigenMatrixColMajor, 0, Eigen::OuterStride<>>;
        StridedMatrixMap active_matrix_map(A.data(), n_rows_filled, n_cols, Eigen::OuterStride<>(total_rows));
        Eigen::Map<EigenVector> b_eigen(b.data(), total_rows);
        EigenVector active_target_vector = b_eigen.head(n_rows_filled);

        EigenVector householder_coeffs(n_cols);

        // In-place Householder QR directly on the active top rows of workspace (0 extra matrix allocations!)
        Eigen::internal::householder_qr_inplace_blocked<StridedMatrixMap, EigenVector>::run(
            active_matrix_map, householder_coeffs
        );

        Eigen::HouseholderSequence<StridedMatrixMap, EigenVector> householder_Q(active_matrix_map, householder_coeffs);
        b_eigen.head(n_rows_filled) = householder_Q.transpose() * active_target_vector;

        active_matrix_map.topRows(n_cols).template triangularView<Eigen::StrictlyLower>().setZero();

        n_rows_filled = n_cols;
    }

    void BlockIncrementalQRAccumulator::appendBlock(
        const Matrix<ColumnMajor>& A_block, const std::vector<double>& b_chunk
    ) {
        const size_t block_rows = A_block.nRows();
        assert(A_block.nColumns() == n_cols);
        assert(b_chunk.size() == block_rows);

        if (block_rows == 0) return;

        std::lock_guard lock(mutex);

        const size_t total_rows = n_cols + n_increment_block_rows;
        Eigen::Map<EigenMatrixColMajor> A_eigen(A.data(), total_rows, n_cols);
        Eigen::Map<EigenVector> b_eigen(b.data(), total_rows);

        Eigen::Map<const EigenMatrixColMajor> A_block_eigen(A_block.data(), block_rows, n_cols);
        Eigen::Map<const EigenVector> b_chunk_eigen(b_chunk.data(), block_rows);

        for (size_t block_row = 0; block_row < block_rows; block_row++) {
            A_eigen.row(n_rows_filled) = A_block_eigen.row(block_row);
            b_eigen(n_rows_filled) = b_chunk_eigen(block_row);

            if (++n_rows_filled == total_rows) {
                flushFullBlock();
            }
        }
    }

    std::vector<double> BlockIncrementalQRAccumulator::solve() {
        finalize();

        const size_t total_rows = n_cols + n_increment_block_rows;

        JGAP_LOG_INFO("Streaming QR: Solving triangular system R_accum * c = y_accum ({}x{})", n_cols, n_cols);

        Eigen::Map<const EigenMatrixColMajor> A_eigen(A.data(), total_rows, n_cols);
        Eigen::Map<const EigenVector> b_eigen(b.data(), total_rows);

        // Exact upper-triangular solve R_accum * c = y_accum
        EigenVector solution_coefficients =
            A_eigen.topRows(n_cols).template triangularView<Eigen::Upper>().solve(b_eigen.head(n_cols));

        return std::vector<double>{
            solution_coefficients.data(), solution_coefficients.data() + solution_coefficients.size()
        };
    }

}
