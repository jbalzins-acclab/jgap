#ifndef JGAP_LINALG_HPP
#define JGAP_LINALG_HPP

#include <memory>
#include <optional>
#include <vector>

#include "jgap/core/Matrix.hpp"
#include "jgap/core/Real.hpp"

namespace jgap::linalg {

    /// Solves min ||A*c - b||_2 using unpivoted Householder QR decomposition (ColumnMajor matrix A).
    std::vector<Real> solveLeastSquaresHouseholderQR(Matrix<ColumnMajor>& A, std::vector<Real>& b);

    /// Solves min ||A*c - b||_2 using Least Squares Conjugate Gradient (RowMajor matrix A).
    std::vector<Real> solveLeastSquaresConjugateGradient(
        Matrix<RowMajor>& A,
        std::vector<Real>& b,
        std::optional<int> max_iterations = std::nullopt,
        std::optional<double> tolerance = std::nullopt
    );

    /// Computes Upper Cholesky factor U such that M = U^T * U (from RowMajor matrix M).
    template<MatrixLayout OutputLayout = MatrixLayout::ColumnMajor>
    Matrix<OutputLayout> choleskyDecomposition(Matrix<RowMajor>& matrix_block);

}

#endif

