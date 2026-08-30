#include "jgap/core/linalg/StreamingQRAccumulator.hpp"

#include <Eigen/Dense>
#include <Eigen/QR>
#include <cassert>
#include <vector>

#include "jgap/core/io/log/CurrentLogger.hpp"

namespace jgap::linalg {

    using EigenMatrixCol = Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>;
    using EigenVectorCol = Eigen::Matrix<Real, Eigen::Dynamic, 1>;

    StreamingQRAccumulator::StreamingQRAccumulator(size_t n_cols, size_t max_chunk_rows) :
        n_cols(n_cols),
        max_chunk_rows(std::max<size_t>(1, max_chunk_rows)),
        buffered_rows(0),
        workspace(n_cols + std::max<size_t>(1, max_chunk_rows), n_cols),
        b_workspace(n_cols + std::max<size_t>(1, max_chunk_rows), 0.0_r) {}

    size_t StreamingQRAccumulator::nCols() const { return n_cols; }

    void StreamingQRAccumulator::flushFullBlock() {
        assert(buffered_rows == max_chunk_rows);
        const size_t total_rows = n_cols + max_chunk_rows;
        Eigen::Map<EigenMatrixCol> ws_map(workspace.flatData().data(), total_rows, n_cols);
        Eigen::Map<EigenVectorCol> b_ws_map(b_workspace.data(), total_rows);

        Eigen::HouseholderQR<EigenMatrixCol> qr(ws_map);
        EigenVectorCol Qt_b = qr.householderQ().transpose() * b_ws_map;

        ws_map.topRows(n_cols) = qr.matrixQR().topRows(n_cols).template triangularView<Eigen::Upper>();
        b_ws_map.head(n_cols) = Qt_b.head(n_cols);
        buffered_rows = 0;
    }

    void StreamingQRAccumulator::finalize() {
        if (buffered_rows == 0) {
            return;
        }
        const size_t total_rows = n_cols + max_chunk_rows;
        const size_t active_rows = n_cols + buffered_rows;
        Eigen::Map<const EigenMatrixCol> ws_map(workspace.flatData().data(), total_rows, n_cols);
        Eigen::Map<const EigenVectorCol> b_ws_map(b_workspace.data(), total_rows);

        EigenMatrixCol active_mat(active_rows, n_cols);
        active_mat.topRows(n_cols) = ws_map.topRows(n_cols);
        active_mat.bottomRows(buffered_rows) = ws_map.block(n_cols, 0, buffered_rows, n_cols);

        EigenVectorCol active_b(active_rows);
        active_b.head(n_cols) = b_ws_map.head(n_cols);
        active_b.tail(buffered_rows) = b_ws_map.segment(n_cols, buffered_rows);

        Eigen::HouseholderQR<EigenMatrixCol> qr(active_mat);
        EigenVectorCol Qt_b = qr.householderQ().transpose() * active_b;

        Eigen::Map<EigenMatrixCol> ws_out(workspace.flatData().data(), total_rows, n_cols);
        Eigen::Map<EigenVectorCol> b_out(b_workspace.data(), total_rows);

        ws_out.topRows(n_cols) = qr.matrixQR().topRows(n_cols).template triangularView<Eigen::Upper>();
        b_out.head(n_cols) = Qt_b.head(n_cols);
        buffered_rows = 0;
    }

    void StreamingQRAccumulator::appendBlock(const Matrix<ColumnMajor>& A_chunk, const std::vector<Real>& b_chunk) {
        const size_t chunk_rows = A_chunk.nRows();
        const size_t C = n_cols;
        assert(A_chunk.nColumns() == C);
        assert(b_chunk.size() == chunk_rows);

        if (chunk_rows == 0) return;

        const size_t total_rows = n_cols + max_chunk_rows;
        Eigen::Map<EigenMatrixCol> ws_map(workspace.flatData().data(), total_rows, C);
        Eigen::Map<EigenVectorCol> b_ws_map(b_workspace.data(), total_rows);

        Eigen::Map<const EigenMatrixCol> A_chunk_map(A_chunk.flatData().data(), chunk_rows, C);
        Eigen::Map<const EigenVectorCol> b_chunk_map(b_chunk.data(), chunk_rows);

        size_t chunk_cursor = 0;
        while (chunk_cursor < chunk_rows) {
            const size_t space = max_chunk_rows - buffered_rows;
            const size_t to_copy = std::min(space, chunk_rows - chunk_cursor);

            ws_map.block(n_cols + buffered_rows, 0, to_copy, C) = A_chunk_map.block(chunk_cursor, 0, to_copy, C);
            b_ws_map.segment(n_cols + buffered_rows, to_copy) = b_chunk_map.segment(chunk_cursor, to_copy);

            buffered_rows += to_copy;
            chunk_cursor += to_copy;

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
        Eigen::Map<const EigenMatrixCol> ws_map(workspace.flatData().data(), total_rows, C);
        Eigen::Map<const EigenVectorCol> b_ws_map(b_workspace.data(), total_rows);

        EigenVectorCol c = ws_map.topRows(C).template triangularView<Eigen::Upper>().solve(b_ws_map.head(C));
        return std::vector<Real>{c.data(), c.data() + c.size()};
    }

}
