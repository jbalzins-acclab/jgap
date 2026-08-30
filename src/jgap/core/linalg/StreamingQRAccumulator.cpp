#include "jgap/core/linalg/StreamingQRAccumulator.hpp"

#include <Eigen/Dense>
#include <Eigen/QR>
#include <cassert>
#include <vector>

#include "jgap/core/io/log/CurrentLogger.hpp"

namespace jgap::linalg {

    // Dynamic-size column-major matrix matching jgap::Matrix<ColumnMajor> layout
    using EigenMatrixColMajor = Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>;
    // Dynamic-size column vector
    using EigenVector = Eigen::Matrix<Real, Eigen::Dynamic, 1>;

    size_t StreamingQRAccumulator::calculateMaxChunkRows(const size_t n_cols, const double approx_ram_limit_gb) {
        const double max_bytes = approx_ram_limit_gb * 1024.0 * 1024.0 * 1024.0;
        const double bytes_per_row = static_cast<double>(n_cols) * sizeof(Real);
        const double mm_bytes = static_cast<double>(n_cols) * bytes_per_row;

        // In-place peak streaming memory breakdown:
        // 1. workspace: (n_cols + B) * n_cols * 8 bytes (sizeof(Real))
        // 2. A_chunk input buffer: B * n_cols * 8 bytes
        // Total Peak RAM = (n_cols + B) * n_cols * 8 bytes + B * n_cols * 8 bytes
        //                = mm_bytes + 2 * B * bytes_per_row
        const double min_required_bytes = mm_bytes + 2.0 * bytes_per_row;
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
        const size_t max_b = static_cast<size_t>(remaining_bytes / (2.0 * bytes_per_row));
        return std::max<size_t>(1, max_b);
    }

    StreamingQRAccumulator::StreamingQRAccumulator(const size_t n_cols, const double approx_ram_limit_gb) :
        n_cols(n_cols),
        max_chunk_rows(calculateMaxChunkRows(n_cols, approx_ram_limit_gb)),
        buffered_rows(0),
        workspace(n_cols + max_chunk_rows, n_cols),
        b_workspace(n_cols + max_chunk_rows, 0.0_r) {}

    size_t StreamingQRAccumulator::nCols() const { return n_cols; }

    void StreamingQRAccumulator::flushFullBlock() {
        assert(buffered_rows == max_chunk_rows);
        const size_t total_rows = n_cols + max_chunk_rows;

        // Map workspace memory without dynamic heap allocation
        Eigen::Map<EigenMatrixColMajor> workspace_matrix_map(workspace.flatData().data(), total_rows, n_cols);
        Eigen::Map<EigenVector> workspace_target_vector_map(b_workspace.data(), total_rows);

        // Vector of Householder scalar coefficients (length = n_cols, only a few KB)
        EigenVector hCoeffs(n_cols);

        // In-place blocked Householder QR directly inside workspace_matrix_map (0 extra matrix allocations!)
        Eigen::internal::householder_qr_inplace_blocked<Eigen::Map<EigenMatrixColMajor>, EigenVector>::run(
            workspace_matrix_map, hCoeffs
        );

        // Apply Q^T in-place to the target vector using the implicit Householder sequence
        Eigen::HouseholderSequence<Eigen::Map<EigenMatrixColMajor>, EigenVector> householder_Q(
            workspace_matrix_map, hCoeffs
        );
        workspace_target_vector_map = householder_Q.transpose() * workspace_target_vector_map;

        // Zero out the strictly lower-triangular part of the top n_cols x n_cols block of R
        workspace_matrix_map.topRows(n_cols).template triangularView<Eigen::StrictlyLower>().setZero();

        buffered_rows = 0;
    }

    void StreamingQRAccumulator::finalize() {
        if (buffered_rows == 0) {
            return;
        }
        const size_t total_rows = n_cols + max_chunk_rows;
        const size_t active_rows = n_cols + buffered_rows;

        // Use explicit OuterStride for the active_rows x n_cols submatrix to guarantee
        // optimal vectorization and avoid complex nested template types in Householder QR.
        using StridedMatrixMap = Eigen::Map<EigenMatrixColMajor, 0, Eigen::OuterStride<>>;
        StridedMatrixMap active_matrix_map(
            workspace.flatData().data(), active_rows, n_cols, Eigen::OuterStride<>(total_rows)
        );
        Eigen::Map<EigenVector> workspace_target_vector_map(b_workspace.data(), total_rows);
        EigenVector active_target_vector = workspace_target_vector_map.head(active_rows);

        EigenVector hCoeffs(n_cols);

        // In-place Householder QR directly on the active top rows of workspace (0 extra matrix allocations!)
        Eigen::internal::householder_qr_inplace_blocked<StridedMatrixMap, EigenVector>::run(active_matrix_map, hCoeffs);

        Eigen::HouseholderSequence<StridedMatrixMap, EigenVector> householder_Q(active_matrix_map, hCoeffs);
        workspace_target_vector_map.head(active_rows) = householder_Q.transpose() * active_target_vector;

        active_matrix_map.topRows(n_cols).template triangularView<Eigen::StrictlyLower>().setZero();

        buffered_rows = 0;
    }

    void StreamingQRAccumulator::appendBlock(const Matrix<ColumnMajor>& A_chunk, const std::vector<Real>& b_chunk) {
        const size_t chunk_rows = A_chunk.nRows();
        const size_t C = n_cols;
        assert(A_chunk.nColumns() == C);
        assert(b_chunk.size() == chunk_rows);

        if (chunk_rows == 0) return;

        const size_t total_rows = n_cols + max_chunk_rows;
        Eigen::Map<EigenMatrixColMajor> workspace_matrix_map(workspace.flatData().data(), total_rows, C);
        Eigen::Map<EigenVector> workspace_target_vector_map(b_workspace.data(), total_rows);

        Eigen::Map<const EigenMatrixColMajor> chunk_matrix_map(A_chunk.flatData().data(), chunk_rows, C);
        Eigen::Map<const EigenVector> chunk_target_vector_map(b_chunk.data(), chunk_rows);

        size_t chunk_cursor = 0;
        while (chunk_cursor < chunk_rows) {
            const size_t space = max_chunk_rows - buffered_rows;
            const size_t to_copy = std::min(space, chunk_rows - chunk_cursor);

            // Copy slice of A_chunk into workspace starting below current R and buffered rows:
            // workspace[n_cols + buffered_rows : n_cols + buffered_rows + to_copy, 0 : C]
            workspace_matrix_map.block(n_cols + buffered_rows, 0, to_copy, C) =
                chunk_matrix_map.block(chunk_cursor, 0, to_copy, C);

            // Copy slice of b_chunk into b_workspace segment
            workspace_target_vector_map.segment(n_cols + buffered_rows, to_copy) =
                chunk_target_vector_map.segment(chunk_cursor, to_copy);

            buffered_rows += to_copy;
            chunk_cursor += to_copy;

            // When buffer is full, perform QR decomposition to absorb chunk into R
            if (buffered_rows == max_chunk_rows) {
                flushFullBlock();
            }
        }
    }

    std::vector<Real> StreamingQRAccumulator::solve() {
        finalize();
        const size_t C = n_cols;
        const size_t total_rows = n_cols + max_chunk_rows;
        JGAP_LOG_INFO("Streaming QR: Solving triangular system R_accum * c = y_accum ({}x{})", C, C);

        Eigen::Map<const EigenMatrixColMajor> workspace_matrix_map(workspace.flatData().data(), total_rows, C);
        Eigen::Map<const EigenVector> workspace_target_vector_map(b_workspace.data(), total_rows);

        // Solve upper-triangular system R_accum * c = y_accum using backward substitution.
        // triangularView<Eigen::Upper>().solve() leverages optimized BLAS/SIMD triangular solvers.
        EigenVector solution_coefficients =
            workspace_matrix_map.topRows(C).template triangularView<Eigen::Upper>().solve(
                workspace_target_vector_map.head(C)
            );

        return std::vector<Real>{
            solution_coefficients.data(), solution_coefficients.data() + solution_coefficients.size()
        };
    }

}
