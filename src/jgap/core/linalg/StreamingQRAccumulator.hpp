#ifndef JGAP_STREAMINGQRACCUMULATOR_HPP
#define JGAP_STREAMINGQRACCUMULATOR_HPP

#include <optional>
#include <vector>

#include "jgap/core/Matrix.hpp"
#include "jgap/core/Real.hpp"

namespace jgap::linalg {

    /// Out-of-core / chunked QR accumulator that pre-allocates a single contiguous (M + max_chunk_rows) x M
    /// workspace and buffers chunks to perform standard Eigen HouseholderQR factorizations in-place.
    class StreamingQRAccumulator {
    public:
        StreamingQRAccumulator(size_t n_cols, double approx_ram_limit_gb);
        ~StreamingQRAccumulator() = default;

        StreamingQRAccumulator(StreamingQRAccumulator&&) noexcept = default;
        StreamingQRAccumulator& operator=(StreamingQRAccumulator&&) noexcept = default;

        StreamingQRAccumulator(const StreamingQRAccumulator&) = default;
        StreamingQRAccumulator& operator=(const StreamingQRAccumulator&) = default;

        Real& operator()(size_t i, size_t j) { return workspace(i, j); }
        const Real& operator()(size_t i, size_t j) const { return workspace(i, j); }

        Matrix<ColumnMajor>& initialWorkspace() { return workspace; }
        const Matrix<ColumnMajor>& initialWorkspace() const { return workspace; }

        std::vector<Real>& initialB() { return b_workspace; }
        const std::vector<Real>& initialB() const { return b_workspace; }

        void appendBlock(const Matrix<ColumnMajor>& A_chunk, const std::vector<Real>& b_chunk);
        std::vector<Real> solve();
        size_t nCols() const;
        size_t targetChunkRows() const { return max_chunk_rows; }

        void flush() { finalize(); }

    private:
        static size_t calculateMaxChunkRows(size_t n_cols, double approx_ram_limit_gb);
        void flushFullBlock();
        void finalize();

        size_t n_cols{0};
        size_t max_chunk_rows{1000};
        size_t buffered_rows{0};
        Matrix<ColumnMajor> workspace;
        std::vector<Real> b_workspace;
    };

}

#endif
