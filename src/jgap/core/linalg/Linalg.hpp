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

    /// Out-of-core / chunked QR accumulator that incrementally maintains upper-triangular R and Q^T * b.
    class StreamingQRAccumulator {
    public:
        explicit StreamingQRAccumulator(size_t n_cols);
        ~StreamingQRAccumulator();

        StreamingQRAccumulator(StreamingQRAccumulator&&) noexcept;
        StreamingQRAccumulator& operator=(StreamingQRAccumulator&&) noexcept;

        StreamingQRAccumulator(const StreamingQRAccumulator&) = delete;
        StreamingQRAccumulator& operator=(const StreamingQRAccumulator&) = delete;

        void appendBlock(const Matrix<ColumnMajor>& A_chunk, const std::vector<Real>& b_chunk);
        std::vector<Real> solve() const;
        size_t nCols() const;

    private:
        struct Impl;
        std::unique_ptr<Impl> pimpl;
    };

} // namespace jgap::linalg

#endif // JGAP_LINALG_HPP
