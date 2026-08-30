#include "jgap/experimental/fit/gap/StreamingQrGapFit.hpp"

#include "../../../core/io/log/CurrentLogger.hpp"
#include "jgap/core/UnseqFor.hpp"
#include "jgap/core/linalg/Linalg.hpp"
#include "jgap/core/linalg/StreamingQRAccumulator.hpp"

namespace jgap {
    StreamingQrGapFit::StreamingQrGapFit(const Real jitter, const double approx_ram_limit_gb) :
        QRGapFit(jitter), approx_ram_limit_gb(approx_ram_limit_gb) {}

    std::vector<Real> StreamingQrGapFit::findCoefficients(
        std::vector<ValuePtr<GapComponent>>& gap_components,
        const std::vector<Atoms>& training_data,
        std::vector<EnergyData>& energies_without_external,
        std::vector<Regularization>& sigmas_inverse
    ) {
        size_t C = 0;
        for (const auto& comp: gap_components) {
            C += comp->nSparsePoints();
        }

        linalg::StreamingQRAccumulator accumulator(C, approx_ram_limit_gb);
        const size_t effective_chunk_rows = accumulator.targetChunkRows();

        size_t start_col = 0;
        for (const auto& comp: gap_components) {
            fillU_mm(start_col, start_col, comp, accumulator.initialWorkspace());
            start_col += comp->nSparsePoints();
        }

        struct FrameMeta {
            size_t frame_idx;
            size_t rows;
        };

        std::vector<FrameMeta> frames_meta;
        size_t N = 0;
        for (size_t i = 0; i < training_data.size(); ++i) {
            size_t r = 0;
            if (energies_without_external[i].energy.has_value()) r += 1;
            if (energies_without_external[i].forces.has_value()) r += 3 * training_data[i].nAtoms();
            if (energies_without_external[i].virials.has_value()) r += 6;
            frames_meta.push_back({i, r});
            N += r;
        }

        const double bytes_per_row = static_cast<double>(C) * sizeof(Real);
        const double mm_mb = (static_cast<double>(C) * static_cast<double>(C) * sizeof(Real)) / (1024.0 * 1024.0);
        const double full_nm_mb = (static_cast<double>(N + C) * bytes_per_row) / (1024.0 * 1024.0);
        const double chunk_mb = (static_cast<double>(effective_chunk_rows) * bytes_per_row) / (1024.0 * 1024.0);
        const double est_peak_gb =
            (static_cast<double>(C + effective_chunk_rows) * bytes_per_row + static_cast<double>(effective_chunk_rows) * bytes_per_row) /
            (1024.0 * 1024.0 * 1024.0);

        JGAP_LOG_INFO(
            "Starting Streaming QR fit: columns C = {}, rows N = {}, RAM/row = {:.2f} KB, MxM RAM = {:.2f} MB, full "
            "(N+M)xM RAM = {:.2f} MB, limit = {:.2f} GB, chunk capacity = {} rows ({:.2f} MB), est. peak = {:.2f} GB",
            C,
            N,
            bytes_per_row / 1024.0,
            mm_mb,
            full_nm_mb,
            approx_ram_limit_gb,
            effective_chunk_rows,
            chunk_mb,
            est_peak_gb
        );

        size_t frame_cursor = 0;
        size_t chunk_count = 0;

        while (frame_cursor < frames_meta.size()) {
            size_t chunk_start_frame = frame_cursor;
            size_t chunk_rows = 0;

            while (frame_cursor < frames_meta.size() && (chunk_rows < effective_chunk_rows || chunk_rows == 0)) {
                chunk_rows += frames_meta[frame_cursor].rows;
                frame_cursor++;
            }
            size_t chunk_end_frame = frame_cursor;
            chunk_count++;

            JGAP_LOG_INFO(
                "Streaming QR: Chunk {} (frames {}..{}, rows {} x cols {})",
                chunk_count,
                chunk_start_frame,
                chunk_end_frame - 1,
                chunk_rows,
                C
            );

            Matrix<ColumnMajor> A_chunk(chunk_rows, C);
            std::vector<Real> b_chunk;
            b_chunk.reserve(chunk_rows);

            std::vector<std::pair<size_t, size_t>> chunk_frame_layout;
            size_t cur_row = 0;
            for (size_t f = chunk_start_frame; f < chunk_end_frame; ++f) {
                size_t f_idx = frames_meta[f].frame_idx;
                chunk_frame_layout.push_back({f_idx, cur_row});

                const auto& ed = energies_without_external[f_idx];
                const auto& sig = sigmas_inverse[f_idx];
                if (ed.energy.has_value()) {
                    b_chunk.push_back(ed.energy.value() * sig.energy.value());
                }
                if (ed.forces.has_value()) {
                    for (size_t j = 0; j < ed.forces->size(); ++j) {
                        b_chunk.push_back(ed.forces->at(j).x * sig.forces->at(j).x);
                        b_chunk.push_back(ed.forces->at(j).y * sig.forces->at(j).y);
                        b_chunk.push_back(ed.forces->at(j).z * sig.forces->at(j).z);
                    }
                }
                if (ed.virials.has_value()) {
                    b_chunk.push_back(ed.virials->xx * sig.virials->xx);
                    b_chunk.push_back(ed.virials->xy * sig.virials->xy);
                    b_chunk.push_back(ed.virials->xz * sig.virials->xz);
                    b_chunk.push_back(ed.virials->yy * sig.virials->yy);
                    b_chunk.push_back(ed.virials->yz * sig.virials->yz);
                    b_chunk.push_back(ed.virials->zz * sig.virials->zz);
                }

                cur_row += frames_meta[f].rows;
            }

            unseqForEach(chunk_frame_layout.begin(), chunk_frame_layout.end(), [&](const auto& layout) {
                size_t f_idx = layout.first;
                size_t start_r = layout.second;
                fillInverseSigmaLK_NM(
                    gap_components,
                    training_data[f_idx],
                    energies_without_external[f_idx],
                    sigmas_inverse[f_idx],
                    A_chunk,
                    start_r
                );
            });

            accumulator.appendBlock(A_chunk, b_chunk);
        }

        return accumulator.solve();
    }
}
