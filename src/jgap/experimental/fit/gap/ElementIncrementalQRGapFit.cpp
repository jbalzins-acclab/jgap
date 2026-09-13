#include "jgap/experimental/fit/gap/ElementIncrementalQRGapFit.hpp"

#include <algorithm>
#include <memory>
#include <set>
#include <vector>

#include "../../../core/io/log/CurrentLogger.hpp"
#include "jgap/core/UnseqFor.hpp"
#include "jgap/core/linalg/BlockIncrementalQRAccumulator.hpp"
#include "jgap/core/linalg/Linalg.hpp"

namespace jgap {

    ElementIncrementalQRGapFit::ElementIncrementalQRGapFit(const double jitter, const double approx_ram_limit_gb) :
        QRGapFit(jitter), approx_ram_limit_gb(approx_ram_limit_gb) {}

    std::vector<double> ElementIncrementalQRGapFit::findCoefficients(
        std::vector<ValuePtr<GapComponent>>& gap_components,
        const std::vector<Atoms>& training_data,
        std::vector<EnergyData>& energies_without_external,
        std::vector<Regularization>& sigmas_inverse
    ) {
        size_t n_col = 0;
        for (const auto& comp: gap_components) {
            n_col += comp->nSparsePoints();
        }

        const auto matrix_memory_space = linalg::BlockIncrementalQRAccumulator::allocateMatrixMemory(
            n_col,
            approx_ram_limit_gb
        );
        const size_t total_rows = matrix_memory_space.size() / n_col;
        SharedArray<double> b(total_rows);

        std::vector<ValuePtr<GapComponent>> unused_gap_components = std::move(gap_components);

        std::set<Species> all_species;
        std::vector<std::set<Species>> species_per_entry(training_data.size());
        for (size_t i = 0; i < training_data.size(); i++) {
            for (size_t j = 0; j < training_data[i].nAtoms(); j++) {
                species_per_entry[i].insert(training_data[i].getSpecies()[j]);
                all_species.insert(training_data[i].getSpecies()[j]);
            }
        }
        for (const auto& comp: unused_gap_components) {
            for (const auto& sp: comp->nonZeroCovarianceFor()) {
                all_species.insert(sp);
            }
        }

        const double bytes_per_row = static_cast<double>(n_col) * sizeof(double);
        const double workspace_mb = (static_cast<double>(total_rows) * bytes_per_row) / (1024.0 * 1024.0);

        std::string species_str;
        for (const auto& sp: all_species) {
            if (!species_str.empty()) species_str += ", ";
            species_str += sp.symbol();
        }

        JGAP_LOG_INFO(
            "Starting Element-Incremental QR fit: total columns M = {}, workspace rows = {} ({:.2f} MB), "
            "RAM limit = {:.2f} GB, species ({}) = [{}], entries = {}",
            n_col,
            total_rows,
            workspace_mb,
            approx_ram_limit_gb,
            all_species.size(),
            species_str,
            training_data.size()
        );

        std::vector<bool> used_entries(training_data.size(), false);
        size_t total_entries_processed = 0;

        size_t current_M = 0;
        std::vector<ValuePtr<GapComponent>> used_gap_components{};
        std::set<Species> used_species{};
        std::unique_ptr<linalg::BlockIncrementalQRAccumulator> last_accumulator;

        size_t elem_idx = 0;
        for (const Species& current_element: all_species) {
            elem_idx++;
            const bool is_last_element = (elem_idx == all_species.size());
            const std::string elem_name = current_element.symbol();

            JGAP_LOG_INFO(
                "[{}/{}] Processing element '{}' (current active M = {}/{})",
                elem_idx,
                all_species.size(),
                elem_name,
                current_M,
                n_col
            );

            // 0. Decompose single-element components for current_element
            std::vector<ValuePtr<GapComponent>> single_element_components;
            for (int i = 0; i < unused_gap_components.size(); i++) {
                std::set<Species> component_species = unused_gap_components[i]->nonZeroCovarianceFor();

                if (component_species.size() == 1 && component_species.contains(current_element)) {
                    single_element_components.push_back(unused_gap_components[i]);
                    unused_gap_components.erase(unused_gap_components.begin() + i);
                    i--;
                }
            }

            size_t M_single = 0;
            if (!single_element_components.empty()) {
                for (const auto& comp: single_element_components) {
                    M_single += comp->nSparsePoints();
                }

                const size_t single_rows = total_rows - current_M;
                const size_t contagious_MxM_block_elements = current_M * current_M;
                auto free_matrix_memory = matrix_memory_space.subspace(contagious_MxM_block_elements);

                Matrix<ColumnMajor> A_single(free_matrix_memory, single_rows, M_single);
                SharedArray<double> b_single(single_rows);

                size_t single_col_cursor = 0;
                for (const auto& comp: single_element_components) {
                    auto U = U_MM(comp);
                    const size_t n = U.nRows();

                    for (size_t j = 0; j < n; j++) {
                        for (size_t r = 0; r < single_col_cursor; ++r) A_single(r, single_col_cursor + j) = 0.0;
                        for (size_t i = 0; i < n; i++) {
                            A_single(single_col_cursor + i, single_col_cursor + j) = U(i, j);
                        }
                        for (size_t r = single_col_cursor + n; r < A_single.nRows(); ++r) A_single(r, single_col_cursor + j) = 0.0;
                    }

                    single_col_cursor += n;
                }

                std::vector<size_t> single_element_entry_indices;
                for (size_t i = 0; i < species_per_entry.size(); i++) {
                    if (!used_entries[i] && species_per_entry[i].size() == 1 && species_per_entry[i].contains(current_element)) {
                        single_element_entry_indices.push_back(i);
                        used_entries[i] = true;
                    }
                }

                JGAP_LOG_INFO(
                    "  Decomposing {} single-element components for '{}' (M_single = {}) with {} single-element training entries",
                    single_element_components.size(),
                    elem_name,
                    M_single,
                    single_element_entry_indices.size()
                );

                linalg::BlockIncrementalQRAccumulator single_accumulator(A_single, b_single, M_single);

                std::atomic<size_t> single_counter(0);
                const size_t total_single = single_element_entry_indices.size();
                const size_t log_interval_single = std::max(total_single / 10, 1uz);

                unseqForEach(
                    single_element_entry_indices.begin(),
                    single_element_entry_indices.end(),
                    [&](const auto entry_idx) {
                        auto matrix_block = formInverseSigmaLK_NMForEntry(
                            single_element_components,
                            training_data[entry_idx],
                            energies_without_external[entry_idx],
                            sigmas_inverse[entry_idx]
                        );
                        auto target_chunk = formTargetVectorBForEntry(
                            energies_without_external[entry_idx],
                            sigmas_inverse[entry_idx]
                        );

                        single_accumulator.appendBlock(matrix_block, target_chunk);

                        const size_t progress = ++single_counter;
                        if (progress % log_interval_single == 0 || progress == total_single) {
                            JGAP_LOG_INFO(
                                "    Single-element '{}' streaming progress: {} of {} entries ({}%)",
                                elem_name,
                                progress,
                                total_single,
                                progress * 100 / total_single
                            );
                        }
                    }
                );

                single_accumulator.finalize();
                A_single.makeContiguous();

                std::copy_n(b_single.data(), M_single, b.data() + current_M);
                total_entries_processed += single_element_entry_indices.size();
            }

            used_species.insert(current_element);

            // 1. Find cross-components whose species are all in used_species
            std::vector<ValuePtr<GapComponent>> current_cross_components;
            for (int i = 0; i < unused_gap_components.size(); i++) {
                const auto comp_species = unused_gap_components[i]->nonZeroCovarianceFor();
                bool all_in = !comp_species.empty() && std::ranges::all_of(comp_species, [&](const Species& s) {
                    return used_species.contains(s);
                });
                if (all_in) {
                    current_cross_components.push_back(unused_gap_components[i]);
                    unused_gap_components.erase(unused_gap_components.begin() + i);
                    i--;
                }
            }

            size_t M_cross = 0;
            for (const auto& comp: current_cross_components) {
                M_cross += comp->nSparsePoints();
            }

            const size_t M_new = current_M + M_single + M_cross;
            if (M_new == 0) {
                continue;
            }

            if (!current_cross_components.empty()) {
                JGAP_LOG_INFO(
                    "  Found {} cross-components for accumulated species (M_cross = {})",
                    current_cross_components.size(),
                    M_cross
                );
            }

            JGAP_LOG_INFO(
                "  Unstacking active system to cumulative M = {}/{} (prev_M = {}, single_M = {}, cross_M = {})",
                M_new,
                n_col,
                current_M,
                M_single,
                M_cross
            );

            // 2. Unstack matrices to new layout
            Matrix<ColumnMajor> A_active(matrix_memory_space, total_rows, M_new);

            if (M_single > 0) {
                Matrix<ColumnMajor> R_single_src(matrix_memory_space.subspace(current_M * current_M), M_single, M_single);
                Matrix<ColumnMajor>::unstack(R_single_src, A_active, current_M, current_M);
            }

            if (current_M > 0) {
                Matrix<ColumnMajor> R_prev_src(matrix_memory_space, current_M, current_M);
                Matrix<ColumnMajor>::unstack(R_prev_src, A_active, 0, 0);
            }

            for (auto& comp: single_element_components) {
                used_gap_components.push_back(std::move(comp));
            }

            // 3. Form U_MM for all new cross-components and add to current matrix
            size_t cross_col_cursor = current_M + M_single;
            for (const auto& comp: current_cross_components) {
                auto U = U_MM(comp);
                const size_t n = U.nRows();

                for (size_t j = 0; j < n; j++) {
                    for (size_t r = 0; r < cross_col_cursor; ++r) A_active(r, cross_col_cursor + j) = 0.0;
                    for (size_t i = 0; i < n; i++) {
                        A_active(cross_col_cursor + i, cross_col_cursor + j) = U(i, j);
                    }
                    for (size_t r = cross_col_cursor + n; r < A_active.nRows(); ++r) A_active(r, cross_col_cursor + j) = 0.0;
                }

                cross_col_cursor += n;
                used_gap_components.push_back(std::move(comp));
            }

            std::fill_n(b.data() + current_M + M_single, M_cross, 0.0);

            // 4. Determine unused entries in training data whose species are all in used_species
            std::vector<size_t> matched_entry_indices;
            for (size_t i = 0; i < training_data.size(); i++) {
                if (!used_entries[i]) {
                    const auto& sp = species_per_entry[i];
                    if (!sp.empty() && std::ranges::all_of(sp, [&](const Species& s) { return used_species.contains(s); })) {
                        matched_entry_indices.push_back(i);
                        used_entries[i] = true;
                    }
                }
            }

            last_accumulator = std::make_unique<linalg::BlockIncrementalQRAccumulator>(
                A_active, b, M_new
            );

            if (!matched_entry_indices.empty()) {
                JGAP_LOG_INFO(
                    "  Streaming {} matched multi-element training structures into active system",
                    matched_entry_indices.size()
                );

                std::atomic<size_t> matched_counter(0);
                const size_t total_matched = matched_entry_indices.size();
                const size_t log_interval_matched = std::max(total_matched / 10, 1uz);

                unseqForEach(
                    matched_entry_indices.begin(),
                    matched_entry_indices.end(),
                    [&](const auto entry_idx) {
                        auto matrix_block = formInverseSigmaLK_NMForEntry(
                            used_gap_components,
                            training_data[entry_idx],
                            energies_without_external[entry_idx],
                            sigmas_inverse[entry_idx]
                        );
                        auto target_chunk = formTargetVectorBForEntry(
                            energies_without_external[entry_idx],
                            sigmas_inverse[entry_idx]
                        );

                        last_accumulator->appendBlock(matrix_block, target_chunk);

                        const size_t progress = ++matched_counter;
                        if (progress % log_interval_matched == 0 || progress == total_matched) {
                            JGAP_LOG_INFO(
                                "    Matched multi-element streaming progress (up to '{}'): {} of {} entries ({}%)",
                                elem_name,
                                progress,
                                total_matched,
                                progress * 100 / total_matched
                            );
                        }
                    }
                );

                total_entries_processed += matched_entry_indices.size();
            }

            last_accumulator->finalize();

            JGAP_LOG_INFO(
                "  Element '{}' finalized: cumulative M = {}/{}, dataset progress = {} of {} entries ({}%)",
                elem_name,
                M_new,
                n_col,
                total_entries_processed,
                training_data.size(),
                training_data.empty() ? 100 : (total_entries_processed * 100 / training_data.size())
            );

            if (!is_last_element) {
                A_active.makeContiguous();
            }

            current_M = M_new;
        }

        if (!last_accumulator) {
            JGAP_LOG_AND_THROW("No components or data processed in ElementIncrementalQRGapFit.");
        }

        JGAP_LOG_INFO(
            "Element-Incremental QR fit complete across all {} species. Solving final system ({}x{})...",
            all_species.size(),
            current_M,
            current_M
        );

        gap_components = std::move(used_gap_components);
        return last_accumulator->solve();
    }
}
