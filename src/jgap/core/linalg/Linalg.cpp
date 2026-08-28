#include "Linalg.hpp"

#include <Eigen/Dense>
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/QR>
#include <cassert>
#include <cmath>

#include "jgap/core/io/log/CurrentLogger.hpp"

namespace jgap::linalg {

    using EigenMatrixCol = Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>;
    using EigenMatrixRow = Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
    using EigenVectorCol = Eigen::Matrix<Real, Eigen::Dynamic, 1>;

    std::vector<Real> solveLeastSquaresHouseholderQR(Matrix<ColumnMajor>& A, std::vector<Real>& b) {
        Eigen::Map<EigenMatrixCol> A_map(A.flatData().data(), A.nRows(), A.nColumns());
        Eigen::Map<EigenVectorCol> b_map(b.data(), b.size());

        JGAP_LOG_DEBUG("Init Eigen::HouseholderQR");
        const Eigen::HouseholderQR<Eigen::Ref<EigenMatrixCol>> qr(A_map);

        JGAP_LOG_DEBUG("Q^t");
        auto Qt = qr.householderQ().transpose();

        JGAP_LOG_DEBUG("Q^t * b");
        auto Qt_b = Qt * b_map;

        JGAP_LOG_DEBUG("R");
        auto R = qr.matrixQR().topLeftCorner(A.nColumns(), A.nColumns());

        JGAP_LOG_DEBUG("R^-1 * Q^t_b");
        EigenVectorCol c = R.triangularView<Eigen::Upper>().solve(Qt_b.head(A.nColumns()));

        Real b_norm = b_map.norm();
        if (b_norm > 0.0_r && A.nRows() > A.nColumns()) {
            Real rel_err = Qt_b.tail(A.nRows() - A.nColumns()).norm() / b_norm;
            JGAP_LOG_INFO("QR solver finished: relative error = {}", rel_err);
        }

        return std::vector<Real>{c.data(), c.data() + c.size()};
    }

    std::vector<Real> solveLeastSquaresConjugateGradient(
        Matrix<RowMajor>& A, std::vector<Real>& b, std::optional<int> max_iterations, std::optional<double> tolerance
    ) {
        const size_t R = A.nRows();
        const size_t C = A.nColumns();

        Eigen::Map<EigenMatrixRow> A_map(A.flatData().data(), R, C);
        Eigen::Map<EigenVectorCol> b_map(b.data(), b.size());

        JGAP_LOG_INFO("Solving least squares via LeastSquaresConjugateGradient directly on A ({}x{})", R, C);

        Eigen::LeastSquaresConjugateGradient<EigenMatrixRow> lscg;
        if (max_iterations.has_value()) {
            lscg.setMaxIterations(max_iterations.value());
        }
        if (tolerance.has_value()) {
            lscg.setTolerance(tolerance.value());
        }

        lscg.compute(A_map);
        if (lscg.info() != Eigen::Success) {
            JGAP_LOG_ERROR("LeastSquaresConjugateGradient setup failed", true);
        }

        EigenVectorCol c = lscg.solve(b_map);

        JGAP_LOG_INFO("LSCG finished: iterations = {}, relative error = {}", lscg.iterations(), lscg.error());

        return std::vector<Real>{c.data(), c.data() + c.size()};
    }

    template<>
    Matrix<ColumnMajor> choleskyDecomposition<MatrixLayout::ColumnMajor>(Matrix<RowMajor>& matrix_block) {
        Eigen::Map<EigenMatrixRow> map(matrix_block.flatData().data(), matrix_block.nRows(), matrix_block.nColumns());

        Eigen::LLT<Eigen::Map<EigenMatrixRow>> llt(map);
        if (llt.info() != Eigen::Success) {
            JGAP_LOG_ERROR("Cholesky decomposition failed: matrix not positive definite", true);
        }

        Matrix<ColumnMajor> result(matrix_block.nRows(), matrix_block.nColumns());
        Eigen::Map<EigenMatrixCol> res_map(result.flatData().data(), result.nRows(), result.nColumns());
        res_map = llt.matrixU();
        return result;
    }

    template<>
    Matrix<RowMajor> choleskyDecomposition<MatrixLayout::RowMajor>(Matrix<RowMajor>& matrix_block) {
        Eigen::Map<EigenMatrixRow> map(matrix_block.flatData().data(), matrix_block.nRows(), matrix_block.nColumns());

        Eigen::LLT<Eigen::Map<EigenMatrixRow>> llt(map);
        if (llt.info() != Eigen::Success) {
            JGAP_LOG_ERROR("Cholesky decomposition failed: matrix not positive definite", true);
        }

        Matrix<RowMajor> result(matrix_block.nRows(), matrix_block.nColumns());
        Eigen::Map<EigenMatrixRow> res_map(result.flatData().data(), result.nRows(), result.nColumns());
        res_map = llt.matrixU();
        return result;
    }

    struct StreamingQRAccumulator::Impl {
        size_t n_cols;
        EigenMatrixCol R_accum;
        EigenVectorCol y_accum;

        explicit Impl(size_t n_cols) :
            n_cols(n_cols), R_accum(EigenMatrixCol::Zero(n_cols, n_cols)), y_accum(EigenVectorCol::Zero(n_cols)) {}
    };

    StreamingQRAccumulator::StreamingQRAccumulator(size_t n_cols) : pimpl(std::make_unique<Impl>(n_cols)) {}

    StreamingQRAccumulator::~StreamingQRAccumulator() = default;
    StreamingQRAccumulator::StreamingQRAccumulator(StreamingQRAccumulator&&) noexcept = default;
    StreamingQRAccumulator& StreamingQRAccumulator::operator=(StreamingQRAccumulator&&) noexcept = default;

    size_t StreamingQRAccumulator::nCols() const { return pimpl->n_cols; }

    void StreamingQRAccumulator::appendBlock(const Matrix<ColumnMajor>& A_chunk, const std::vector<Real>& b_chunk) {
        const size_t chunk_rows = A_chunk.nRows();
        const size_t C = pimpl->n_cols;
        assert(A_chunk.nColumns() == C);
        assert(b_chunk.size() == chunk_rows);

        // A_chunk and b_chunk mapping
        Eigen::Map<const EigenMatrixCol> A_chunk_map(A_chunk.flatData().data(), chunk_rows, C);
        Eigen::Map<const EigenVectorCol> b_chunk_map(b_chunk.data(), chunk_rows);

        EigenMatrixCol M_stacked(C + chunk_rows, C);
        M_stacked.topRows(C) = pimpl->R_accum;
        M_stacked.bottomRows(chunk_rows) = A_chunk_map;

        EigenVectorCol b_stacked(C + chunk_rows);
        b_stacked.head(C) = pimpl->y_accum;
        b_stacked.tail(chunk_rows) = b_chunk_map;

        Eigen::HouseholderQR<EigenMatrixCol> qr(M_stacked);
        EigenVectorCol Qt_b = qr.householderQ().transpose() * b_stacked;

        pimpl->R_accum = qr.matrixQR().topRows(C).template triangularView<Eigen::Upper>();
        pimpl->y_accum = Qt_b.head(C);
    }

    std::vector<Real> StreamingQRAccumulator::solve() const {
        const size_t C = pimpl->n_cols;
        JGAP_LOG_INFO("Streaming QR: Solving triangular system R_accum * c = y_accum ({}x{})", C, C);
        EigenVectorCol c = pimpl->R_accum.template triangularView<Eigen::Upper>().solve(pimpl->y_accum);
        return std::vector<Real>{c.data(), c.data() + c.size()};
    }

}
