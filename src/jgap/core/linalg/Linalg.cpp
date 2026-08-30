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

}
