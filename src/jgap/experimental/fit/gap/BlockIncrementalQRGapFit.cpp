#include "jgap/experimental/fit/gap/BlockIncrementalQRGapFit.hpp"

#include <numeric>

#include "../../../core/io/log/CurrentLogger.hpp"
#include "jgap/core/UnseqFor.hpp"
#include "jgap/core/linalg/BlockIncrementalQRAccumulator.hpp"
#include "jgap/core/linalg/Linalg.hpp"

namespace jgap {
    BlockIncrementalQRGapFit::BlockIncrementalQRGapFit(const double jitter, const double approx_ram_limit_gb) :
        QRGapFit(jitter), approx_ram_limit_gb(approx_ram_limit_gb) {}

    void BlockIncrementalQRGapFit::findCoefficients(
        GapPotential& to_be_fit,
        const std::vector<Atoms>& training_data,
        std::vector<EnergyData>& energies_without_external,
        std::vector<Regularization>& sigmas_inverse
    ) {
        size_t n_cols = 0;
        for (const auto& comp: to_be_fit.components) {
            n_cols += comp->nSparsePoints();
        }

        auto A_memory = linalg::BlockIncrementalQRAccumulator::allocateMatrixMemory(n_cols, approx_ram_limit_gb);
        Matrix<ColumnMajor> A(A_memory, n_cols);
        SharedArray<double> b(A.nRows());

        size_t start_col = 0;
        for (const auto& comp: to_be_fit.components) {
            auto U = U_MM(comp);
            const size_t n = U.nRows();
            for (size_t i = 0; i < n; i++) {
                for (size_t j = 0; j < n; j++) {
                    A(start_col + i, start_col + j) = U(i, j);
                }
            }
            start_col += n;
        }

        linalg::BlockIncrementalQRAccumulator accumulator(A_memory, b, n_cols, n_cols);
        const size_t increment_block_rows = accumulator.nIncrementBlockRows();

        const double bytes_per_row = static_cast<double>(n_cols) * sizeof(double);
        const double workspace_mb =
            (static_cast<double>(n_cols + increment_block_rows) * bytes_per_row) / (1024.0 * 1024.0);
        const double chunk_mb = (static_cast<double>(increment_block_rows) * bytes_per_row) / (1024.0 * 1024.0);

        JGAP_LOG_INFO(
            "Starting Block Incremental QR fit: columns M = {}, increment block = {} rows ({:.2f} MB), "
            "workspace = {:.2f} MB, RAM limit = {:.2f} GB, entries = {}",
            n_cols,
            increment_block_rows,
            chunk_mb,
            workspace_mb,
            approx_ram_limit_gb,
            training_data.size()
        );

        std::vector<size_t> entry_indices(training_data.size());
        std::iota(entry_indices.begin(), entry_indices.end(), 0);

        std::atomic<size_t> counter{0};
        const size_t total_entries = entry_indices.size();
        const size_t log_step = std::max(total_entries / 20, 1uz);

        unseqForEach(
            entry_indices.begin(),
            entry_indices.end(),
            [&](const size_t i) {
                auto A_entry = formInverseSigmaLK_NMForEntry(
                    to_be_fit.components,
                    training_data[i],
                    energies_without_external[i],
                    sigmas_inverse[i]
                );
                auto b_entry = formTargetVectorBForEntry(
                    energies_without_external[i],
                    sigmas_inverse[i]
                );

                accumulator.appendBlock(A_entry, b_entry);

                size_t progress = ++counter;
                if (progress % log_step == 0 || progress == total_entries) {
                    JGAP_LOG_INFO(
                        "Structure accumulation progress: {} of {} ({}%)",
                        progress,
                        total_entries,
                        progress * 100 / total_entries
                    );
                }
            }
        );

        const auto c = accumulator.solve();
        to_be_fit.setCoefficients(c);
    }
}
