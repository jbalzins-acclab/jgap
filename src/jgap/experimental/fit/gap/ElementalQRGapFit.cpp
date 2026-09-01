#include "jgap/experimental/fit/gap/ElementalQRGapFit.hpp"

#include "../../../core/io/log/CurrentLogger.hpp"
#include "jgap/core/UnseqFor.hpp"
#include "jgap/core/linalg/Linalg.hpp"
#include "jgap/core/linalg/StreamingQRAccumulator.hpp"

namespace jgap {

    namespace {
        struct FrameMeta {
            size_t frame_idx;
            size_t rows;
        };
    } // namespace

    ElementalQRGapFit::ElementalQRGapFit(const Real jitter, const double approx_ram_limit_gb) :
        QRGapFit(jitter), approx_ram_limit_gb(approx_ram_limit_gb) {}

    void ElementalQRGapFit::streamFramesIntoAccumulator(
        linalg::StreamingQRAccumulator& accumulator,
        const std::vector<ValuePtr<GapComponent>>& components,
        const std::vector<size_t>& frame_indices,
        const std::vector<Atoms>& training_data,
        const std::vector<EnergyData>& energies_without_external,
        const std::vector<Regularization>& sigmas_inverse,
        const std::string& group_label
    ) const {
        if (frame_indices.empty() || components.empty()) {
            return;
        }

        size_t C = 0;
        for (const auto& comp: components) {
            C += comp->nSparsePoints();
        }

        const size_t effective_chunk_rows = accumulator.targetChunkRows();

        struct LocalFrameMeta {
            size_t frame_idx;
            size_t rows;
        };

        std::vector<LocalFrameMeta> frames_meta;
        size_t N = 0;
        for (size_t f_idx: frame_indices) {
            size_t r = 0;
            if (energies_without_external[f_idx].energy.has_value()) r += 1;
            if (energies_without_external[f_idx].forces.has_value()) r += 3 * training_data[f_idx].nAtoms();
            if (energies_without_external[f_idx].virials.has_value()) r += 6;
            frames_meta.push_back({f_idx, r});
            N += r;
        }

        const double bytes_per_row = static_cast<double>(C) * sizeof(Real);
        const double mm_mb = (static_cast<double>(C) * static_cast<double>(C) * sizeof(Real)) / (1024.0 * 1024.0);
        const double full_nm_mb = (static_cast<double>(N + C) * bytes_per_row) / (1024.0 * 1024.0);
        const double chunk_mb = (static_cast<double>(effective_chunk_rows) * bytes_per_row) / (1024.0 * 1024.0);
        const double est_peak_gb =
            (static_cast<double>(C + effective_chunk_rows) * bytes_per_row +
             static_cast<double>(effective_chunk_rows) * bytes_per_row) /
            (1024.0 * 1024.0 * 1024.0);

        JGAP_LOG_INFO(
            "[{}] Streaming QR fit: columns C = {}, rows N = {}, RAM/row = {:.2f} KB, MxM RAM = {:.2f} MB, full "
            "(N+M)xM RAM = {:.2f} MB, limit = {:.2f} GB, chunk capacity = {} rows ({:.2f} MB), est. peak = {:.2f} GB",
            group_label,
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
            size_t chunk_start = frame_cursor;
            size_t chunk_rows = 0;

            while (frame_cursor < frames_meta.size()) {
                size_t next_rows = frames_meta[frame_cursor].rows;
                if (chunk_rows > 0 && chunk_rows + next_rows > effective_chunk_rows) {
                    break;
                }
                chunk_rows += next_rows;
                frame_cursor++;
            }
            size_t chunk_end = frame_cursor;
            chunk_count++;

            JGAP_LOG_INFO(
                "[{}] Chunk {} (frames {}..{}, rows {} x cols {})",
                group_label,
                chunk_count,
                chunk_start,
                chunk_end - 1,
                chunk_rows,
                C
            );

            Matrix<ColumnMajor> A_chunk(chunk_rows, C);
            std::vector<Real> b_chunk;
            b_chunk.reserve(chunk_rows);

            std::vector<std::pair<size_t, size_t>> chunk_frame_layout;
            size_t cur_row = 0;
            for (size_t f = chunk_start; f < chunk_end; ++f) {
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
                    components,
                    training_data[f_idx],
                    energies_without_external[f_idx],
                    sigmas_inverse[f_idx],
                    A_chunk,
                    start_r
                );
            });

            accumulator.appendBlock(A_chunk, b_chunk);
        }

        accumulator.flush();
    }

    std::vector<Real> ElementalQRGapFit::findCoefficients(
        std::vector<ValuePtr<GapComponent>>& gap_components,
        const std::vector<Atoms>& training_data,
        std::vector<EnergyData>& energies_without_external,
        std::vector<Regularization>& sigmas_inverse
    ) {
        // Collect all unique single species present in the training data
        std::set<Species> unique_species_set;
        for (const auto& atoms: training_data) {
            for (const auto& sp: atoms.getSpecies()) {
                unique_species_set.insert(sp);
            }
        }
        std::vector<Species> species_list(unique_species_set.begin(), unique_species_set.end());

        const size_t K = species_list.size();
        std::vector<std::vector<ValuePtr<GapComponent>>> single_species_components(K);
        std::vector<ValuePtr<GapComponent>> cross_components;

        for (auto& comp: gap_components) {
            std::set<Species> req = comp->nonZeroCovarianceFor();
            std::optional<size_t> matched_idx = std::nullopt;
            if (req.size() == 1) {
                for (size_t k = 0; k < K; ++k) {
                    if (req.contains(species_list[k])) {
                        matched_idx = k;
                        break;
                    }
                }
            }
            if (matched_idx.has_value()) {
                single_species_components[*matched_idx].push_back(comp);
            } else {
                cross_components.push_back(comp);
            }
        }

        // Reorder gap_components so they are grouped as: [G_0, G_1, ..., G_{K-1}, G_cross]
        gap_components.clear();
        for (size_t k = 0; k < K; ++k) {
            for (auto& comp: single_species_components[k]) {
                gap_components.push_back(comp);
            }
        }
        for (auto& comp: cross_components) {
            gap_components.push_back(comp);
        }

        // Classify frames: single-species frames vs cross-species frames
        std::vector<std::vector<size_t>> single_species_frames(K);
        std::vector<size_t> cross_frames;

        for (size_t i = 0; i < training_data.size(); ++i) {
            const auto& spec_vec = training_data[i].getSpecies();
            std::set<Species> frame_species(spec_vec.begin(), spec_vec.end());
            std::optional<size_t> matched_idx = std::nullopt;
            if (frame_species.size() == 1) {
                for (size_t k = 0; k < K; ++k) {
                    if (frame_species.contains(species_list[k])) {
                        matched_idx = k;
                        break;
                    }
                }
            }
            if (matched_idx.has_value()) {
                single_species_frames[*matched_idx].push_back(i);
            } else {
                cross_frames.push_back(i);
            }
        }

        // Assemble full accumulator
        size_t total_C = 0;
        for (const auto& comp: gap_components) {
            total_C += comp->nSparsePoints();
        }

        size_t total_N = 0;
        for (size_t i = 0; i < training_data.size(); ++i) {
            if (energies_without_external[i].energy.has_value()) total_N += 1;
            if (energies_without_external[i].forces.has_value()) total_N += 3 * training_data[i].nAtoms();
            if (energies_without_external[i].virials.has_value()) total_N += 6;
        }

        const double bytes_per_row = static_cast<double>(total_C) * sizeof(Real);
        const double mm_mb =
            (static_cast<double>(total_C) * static_cast<double>(total_C) * sizeof(Real)) / (1024.0 * 1024.0);
        const double full_nm_mb = (static_cast<double>(total_N + total_C) * bytes_per_row) / (1024.0 * 1024.0);

        JGAP_LOG_INFO(
            "Elemental QR: columns C = {}, rows N = {}, RAM/row = {:.2f} KB, MxM RAM = {:.2f} MB, full (N+M)xM RAM = "
            "{:.2f} MB, limit = {:.2f} GB. Classified {} single-species frames and {} cross frames.",
            total_C,
            total_N,
            bytes_per_row / 1024.0,
            mm_mb,
            full_nm_mb,
            approx_ram_limit_gb,
            training_data.size() - cross_frames.size(),
            cross_frames.size()
        );

        linalg::StreamingQRAccumulator full_accum(total_C, approx_ram_limit_gb);

        size_t global_col_offset = 0;
        for (size_t k = 0; k < K; ++k) {
            size_t C_k = 0;
            for (const auto& comp: single_species_components[k]) {
                C_k += comp->nSparsePoints();
            }

            if (C_k > 0) {
                // Scope accum_k so its workspace and chunk buffers are freed immediately
                {
                    linalg::StreamingQRAccumulator accum_k(C_k, approx_ram_limit_gb);
                    size_t start_col = 0;
                    for (const auto& comp: single_species_components[k]) {
                        fillU_mm(start_col, start_col, comp, accum_k.initialWorkspace());
                        start_col += comp->nSparsePoints();
                    }

                    std::string label = "Species " + species_list[k].symbol();
                    streamFramesIntoAccumulator(
                        accum_k,
                        single_species_components[k],
                        single_species_frames[k],
                        training_data,
                        energies_without_external,
                        sigmas_inverse,
                        label
                    );

                    const auto& sub_ws = accum_k.initialWorkspace();
                    const auto& sub_b = accum_k.initialB();
                    for (size_t j = 0; j < C_k; ++j) {
                        for (size_t i = 0; i <= j; ++i) {
                            full_accum(global_col_offset + i, global_col_offset + j) = sub_ws(i, j);
                        }
                    }
                    for (size_t i = 0; i < C_k; ++i) {
                        full_accum.initialB()[global_col_offset + i] = sub_b[i];
                    }
                } // accum_k is destroyed here and all its memory is freed immediately

                single_species_components[k].clear();
                single_species_components[k].shrink_to_fit();
                single_species_frames[k].clear();
                single_species_frames[k].shrink_to_fit();

                global_col_offset += C_k;
            }
        }
        single_species_components.clear();
        single_species_components.shrink_to_fit();
        single_species_frames.clear();
        single_species_frames.shrink_to_fit();

        for (const auto& comp: cross_components) {
            fillU_mm(global_col_offset, global_col_offset, comp, full_accum.initialWorkspace());
            global_col_offset += comp->nSparsePoints();
        }
        cross_components.clear();
        cross_components.shrink_to_fit();

        // Stream cross frames
        if (!cross_frames.empty()) {
            streamFramesIntoAccumulator(
                full_accum,
                gap_components,
                cross_frames,
                training_data,
                energies_without_external,
                sigmas_inverse,
                "Cross"
            );
        }

        return full_accum.solve();
    }

} // namespace jgap
