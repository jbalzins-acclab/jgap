#include "jgap/experimental/fit/gap/ElementIncrementalQRGapFit.hpp"

#include <algorithm>
#include <memory>
#include <set>
#include <vector>

#include "jgap/core/UnseqFor.hpp"
#include "jgap/core/io/log/CurrentLogger.hpp"
#include "jgap/core/linalg/BlockIncrementalQRAccumulator.hpp"
#include "jgap/core/linalg/Linalg.hpp"

namespace jgap {

    ElementIncrementalQRGapFit::ElementIncrementalQRGapFit(const double jitter, const double approx_ram_limit_gb) :
        QRGapFit(jitter), approx_ram_limit_gb(approx_ram_limit_gb) {}

    void ElementIncrementalQRGapFit::findCoefficients(
        GapPotential& to_be_fit,
        const std::vector<Atoms>& training_data,
        std::vector<EnergyData>& energies_without_external,
        std::vector<Regularization>& sigmas_inverse
    ) {
        size_t total_M = 0;
        for (const auto& comp: to_be_fit.components) {
            total_M += comp->nSparsePoints();
        }

        JGAP_LOG_INFO(
            "Starting Element Incremental QR fit: total columns M = {}, components = {}, RAM limit = {:.2f} GB, "
            "entries = {}",
            total_M,
            to_be_fit.components.size(),
            approx_ram_limit_gb,
            training_data.size()
        );

        auto matrix_memory_space = allocateMatrixMemory(to_be_fit.components);

        std::vector<size_t> original_component_indices;
        auto increments = determineIncrements(to_be_fit.components, training_data, original_component_indices);

        std::vector<double> current_b{};
        for (size_t inc_idx = 0; inc_idx < increments.size(); ++inc_idx) {
            IncrementMeta& increment = increments[inc_idx];
            JGAP_LOG_INFO(
                "Starting increment {}/{}: element '{}'", inc_idx + 1, increments.size(), increment.species.symbol()
            );

            size_t current_M = current_b.size();
            size_t contagiousMxM = current_M * current_M;
            auto memory_subspace_for_single_element = matrix_memory_space.subspace(contagiousMxM);

            // ########################## SINGLE ELEMENT ############################
            std::vector<double> b_block{};
            if (!increment.single_element_components.empty()) {
                size_t se_cols = 0;
                for (const auto& comp: increment.single_element_components) {
                    se_cols += comp->nSparsePoints();
                }
                JGAP_LOG_INFO(
                    "Increment {}/{}: fitting {} single-element component(s) (cols = {}, structures = {})",
                    inc_idx + 1,
                    increments.size(),
                    increment.single_element_components.size(),
                    se_cols,
                    increment.single_element_structure_indices.size()
                );

                b_block = singleElementQR(
                    memory_subspace_for_single_element,
                    increment.single_element_components,
                    increment.single_element_structure_indices,
                    training_data,
                    energies_without_external,
                    sigmas_inverse
                );

                for (auto& comp: increment.single_element_components) {
                    to_be_fit.components.push_back(std::move(comp));
                }
                current_b.insert(current_b.end(), b_block.begin(), b_block.end());
            }

            // ########################## NEW CROSS ELEMENTS ############################
            if (!increment.new_multi_element_components.empty()
                || !increment.new_multi_element_structure_indices.empty()) {
                size_t current_n_cols = current_b.size();
                for (const auto& comp: increment.new_multi_element_components) {
                    current_n_cols += comp->nSparsePoints();
                }

                JGAP_LOG_INFO(
                    "Increment {}/{}: restacking matrix and accumulating multi-element data (total cols = {}, new cols "
                    "= {}, structures = {})",
                    inc_idx + 1,
                    increments.size(),
                    current_n_cols,
                    current_n_cols - current_b.size(),
                    increment.new_multi_element_structure_indices.size()
                );

                auto A = restackContagiousBlocks(matrix_memory_space, current_M, b_block.size(), current_n_cols);
                current_M = current_b.size();

                addUMMs(A, current_M, increment.new_multi_element_components);

                auto b = SharedArray<double>(A.nRows());
                for (size_t i = 0; i < current_b.size(); i++) {
                    b[i] = current_b[i];
                }

                for (auto& comp: increment.new_multi_element_components) {
                    to_be_fit.components.push_back(std::move(comp));
                }

                current_b = doQRMakeContagiousAndReturnAccumulatedB(
                    A,
                    b,
                    current_n_cols,
                    to_be_fit.components,
                    increment.new_multi_element_structure_indices,
                    training_data,
                    energies_without_external,
                    sigmas_inverse
                );
            }
        }

        size_t final_M = current_b.size();
        Matrix<ColumnMajor> final_R(matrix_memory_space, final_M, final_M);
        SharedArray<double> final_b(current_b.size());
        for (size_t i = 0; i < current_b.size(); i++) {
            final_b[i] = current_b[i];
        }

        linalg::BlockIncrementalQRAccumulator final_accumulator(final_R, final_b, final_M);
        auto c = final_accumulator.solve();

        to_be_fit.setCoefficients(c);

        std::vector<ValuePtr<GapComponent>> restored_components(to_be_fit.components.size());
        for (size_t k = 0; k < to_be_fit.components.size(); ++k) {
            restored_components[original_component_indices[k]] = std::move(to_be_fit.components[k]);
        }
        to_be_fit.components = std::move(restored_components);

        JGAP_LOG_INFO(
            "Element Incremental QR fit finished: fitted {} coefficients across {} components",
            c.size(),
            to_be_fit.components.size()
        );
    }

    SharedArray<double> ElementIncrementalQRGapFit::allocateMatrixMemory(
        const std::vector<ValuePtr<GapComponent>>& gap_components
    ) const {

        size_t n_col = 0;
        for (const auto& comp: gap_components) {
            n_col += comp->nSparsePoints();
        }

        return linalg::BlockIncrementalQRAccumulator::allocateMatrixMemory(n_col, approx_ram_limit_gb);
    }

    std::vector<ElementIncrementalQRGapFit::IncrementMeta> ElementIncrementalQRGapFit::determineIncrements(
        std::vector<ValuePtr<GapComponent>>& gap_components,
        const std::vector<Atoms>& training_data,
        std::vector<size_t>& original_component_indices
    ) {

        std::set<Species> all_species;
        for (const auto& comp: gap_components) {
            for (const auto& sp: comp->nonZeroCovarianceFor()) {
                all_species.insert(sp);
            }
        }
        std::vector all_species_vec(all_species.begin(), all_species.end());

        std::vector<std::set<Species>> species_per_entry(training_data.size());
        for (size_t i = 0; i < training_data.size(); i++) {
            for (size_t j = 0; j < training_data[i].nAtoms(); j++) {
                species_per_entry[i].insert(training_data[i].getSpecies()[j]);
                all_species.insert(training_data[i].getSpecies()[j]);
            }
        }

        std::vector<bool> entry_used(training_data.size(), false);
        std::vector<bool> component_used(gap_components.size(), false);
        std::vector<IncrementMeta> result(all_species_vec.size());
        std::set<Species> current_species;

        for (int increment = 0; increment < all_species_vec.size(); increment++) {
            Species current_element = all_species_vec[increment];

            current_species.insert(current_element);

            IncrementMeta& meta = result[increment];
            meta.species = current_element;

            // find appropriate training data subsets
            for (size_t i = 0; i < training_data.size(); i++) {
                if (entry_used[i]) continue;
                if (!species_per_entry[i].contains(current_element)) continue;

                bool is_single_species_entry = true;
                bool contains_only_current_species = true;
                for (auto species_in_entry: species_per_entry[i]) {
                    // ignore species not handled by any of the gap components
                    if (!all_species.contains(species_in_entry)) continue;

                    if (species_in_entry != current_element) {
                        is_single_species_entry = false;
                    }
                    if (!current_species.contains(species_in_entry)) {
                        contains_only_current_species = false;
                        break;
                    }
                }

                entry_used[i] = is_single_species_entry || contains_only_current_species;

                if (is_single_species_entry) {
                    meta.single_element_structure_indices.push_back(i);
                } else if (contains_only_current_species) {
                    meta.new_multi_element_structure_indices.push_back(i);
                }
            }

            for (size_t c = 0; c < gap_components.size(); c++) {
                if (component_used[c]) continue;

                auto component_species = gap_components[c]->nonZeroCovarianceFor();

                if (!component_species.contains(current_element)) continue;

                if (component_species.size() == 1) {
                    component_used[c] = true;
                    meta.single_element_components.push_back(std::move(gap_components[c]));
                    original_component_indices.push_back(c);
                }
            }

            for (size_t c = 0; c < gap_components.size(); c++) {
                if (component_used[c]) continue;

                auto component_species = gap_components[c]->nonZeroCovarianceFor();

                if (!component_species.contains(current_element)) continue;

                bool only_current_species = true;
                for (auto cs: component_species) {
                    if (!current_species.contains(cs)) {
                        only_current_species = false;
                        break;
                    }
                }
                if (only_current_species) {
                    component_used[c] = true;
                    meta.new_multi_element_components.push_back(std::move(gap_components[c]));
                    original_component_indices.push_back(c);
                }
            }

            JGAP_LOG_INFO(
                "Determined increment {}/{}: element '{}' -> single-element comps = {}, single structures = {}, "
                "multi-element comps = {}, multi structures = {}",
                increment + 1,
                all_species_vec.size(),
                current_element.symbol(),
                meta.single_element_components.size(),
                meta.single_element_structure_indices.size(),
                meta.new_multi_element_components.size(),
                meta.new_multi_element_structure_indices.size()
            );
        }

        gap_components = {};
        return result;
    }

    std::vector<double> ElementIncrementalQRGapFit::doQRMakeContagiousAndReturnAccumulatedB(
        Matrix<ColumnMajor> A,
        SharedArray<double> b,
        size_t n_cols,
        const std::vector<ValuePtr<GapComponent>>& components,
        const std::vector<size_t>& structure_indices,
        const std::vector<Atoms>& training_data,
        const std::vector<EnergyData>& energy_data,
        const std::vector<Regularization>& sigmas_inverse
    ) {
        linalg::BlockIncrementalQRAccumulator accumulator(A, b, n_cols);

        const size_t increment_block_rows = accumulator.nIncrementBlockRows();
        const double bytes_per_row = static_cast<double>(n_cols) * sizeof(double);
        const double workspace_mb =
            (static_cast<double>(n_cols + increment_block_rows) * bytes_per_row) / (1024.0 * 1024.0);
        const double chunk_mb = (static_cast<double>(increment_block_rows) * bytes_per_row) / (1024.0 * 1024.0);

        JGAP_LOG_INFO(
            "QR Accumulation: columns M = {}, increment block = {} rows ({:.2f} MB), workspace = {:.2f} MB, structures "
            "= {}",
            n_cols,
            increment_block_rows,
            chunk_mb,
            workspace_mb,
            structure_indices.size()
        );

        std::atomic<size_t> counter{0};
        const size_t total_structures = structure_indices.size();
        const size_t log_step = std::max(total_structures / 20, 1uz);

        unseqForEach(structure_indices, [&](size_t i) {
            const auto A_block =
                formInverseSigmaLK_NMForEntry(components, training_data[i], energy_data[i], sigmas_inverse[i]);
            const auto b_block = formTargetVectorBForEntry(energy_data[i], sigmas_inverse[i]);
            accumulator.appendBlock(A_block, b_block);

            size_t progress = ++counter;
            if (progress % log_step == 0 || progress == total_structures) {
                JGAP_LOG_INFO(
                    "Structure accumulation progress: {} of {} ({}%)",
                    progress,
                    total_structures,
                    progress * 100 / total_structures
                );
            }
        });

        accumulator.finalize();

        Matrix<ColumnMajor> contagious_matrix(A.flatData(), n_cols, n_cols);
        for (size_t j = 0; j < n_cols; j++) {
            for (size_t i = 0; i < n_cols; i++) {
                contagious_matrix(i, j) = A(i, j);
            }
        }
        A.flatData().subspace(n_cols * n_cols).fill(0.0);

        std::vector<double> relevant_b_block{};
        relevant_b_block.reserve(n_cols);
        for (size_t i = 0; i < n_cols; i++) {
            relevant_b_block.push_back(b[i]);
        }
        return relevant_b_block;
    }

    std::vector<double> ElementIncrementalQRGapFit::singleElementQR(
        SharedArray<double>& single_element_memory_space,
        const std::vector<ValuePtr<GapComponent>>& components,
        const std::vector<size_t>& single_element_structure_indices,
        const std::vector<Atoms>& training_data,
        const std::vector<EnergyData>& energy_data,
        const std::vector<Regularization>& sigmas_inverse
    ) {
        single_element_memory_space.fill(0.0);

        size_t n_cols = 0;
        for (const auto& comp: components) {
            n_cols += comp->nSparsePoints();
        }

        Matrix<ColumnMajor> A(single_element_memory_space, n_cols);
        SharedArray<double> b(A.nRows());

        size_t current_M{};
        for (size_t comp_idx = 0; comp_idx < components.size(); ++comp_idx) {
            const auto& comp = components[comp_idx];
            JGAP_LOG_INFO(
                "Calculating U_MM for single-element component {}/{} ({} sparse points)",
                comp_idx + 1,
                components.size(),
                comp->nSparsePoints()
            );
            auto U = U_MM(comp);

            const size_t m = U.nRows();
            assert(m == comp->nSparsePoints());

            for (size_t j = 0; j < m; j++) {
                for (size_t i = 0; i < m; i++) {
                    A(current_M + i, current_M + j) = U(i, j);
                }
            }

            current_M += m;
        }

        return doQRMakeContagiousAndReturnAccumulatedB(
            A, b, n_cols, components, single_element_structure_indices, training_data, energy_data, sigmas_inverse
        );
    }

    Matrix<ColumnMajor> ElementIncrementalQRGapFit::restackContagiousBlocks(
        SharedArray<double>& full_memory, size_t former_M, size_t single_element_cols, size_t require_n_cols
    ) {
        JGAP_LOG_INFO(
            "Restacking matrix: former_M = {}, single_element_cols = {}, required_cols = {}",
            former_M,
            single_element_cols,
            require_n_cols
        );
        Matrix<ColumnMajor> restacked(full_memory, require_n_cols);

        // 1. Move single element block backwards in-place
        Matrix<ColumnMajor> contagious_single_element_matrix(
            full_memory.subspace(former_M * former_M), single_element_cols, single_element_cols
        );
        for (size_t j = single_element_cols; j-- > 0;) {
            for (size_t i = single_element_cols; i-- > 0;) {
                restacked(former_M + i, former_M + j) = contagious_single_element_matrix(i, j);
            }
        }

        // 2. Move former species block backwards in-place
        Matrix<ColumnMajor> contagious_former_species_matrix(full_memory, former_M, former_M);
        for (size_t j = former_M; j-- > 0;) {
            for (size_t i = former_M; i-- > 0;) {
                restacked(i, j) = contagious_former_species_matrix(i, j);
            }
        }

        // 3. In-place zero out all off-diagonal and lower rows in restacked
        const size_t total_rows = restacked.nRows();

        // Zero off-diagonal rows in former_M columns
        for (size_t j = 0; j < former_M; ++j) {
            for (size_t i = former_M; i < total_rows; ++i) {
                restacked(i, j) = 0.0;
            }
        }

        // Zero off-diagonal rows in single_element columns
        for (size_t j = 0; j < single_element_cols; ++j) {
            const size_t col = former_M + j;
            for (size_t i = 0; i < former_M; ++i) {
                restacked(i, col) = 0.0;
            }
            for (size_t i = former_M + single_element_cols; i < total_rows; ++i) {
                restacked(i, col) = 0.0;
            }
        }

        // Zero new multi-element columns completely (addUMMs will place U_MM on their diagonal)
        for (size_t col = former_M + single_element_cols; col < require_n_cols; ++col) {
            for (size_t i = 0; i < total_rows; ++i) {
                restacked(i, col) = 0.0;
            }
        }

        return restacked;
    }

    void ElementIncrementalQRGapFit::addUMMs(
        Matrix<ColumnMajor>& A, size_t filled_M, const std::vector<ValuePtr<GapComponent>>& new_multi_element_components
    ) {
        size_t current_M = filled_M;
        for (size_t comp_idx = 0; comp_idx < new_multi_element_components.size(); ++comp_idx) {
            const auto& comp = new_multi_element_components[comp_idx];
            JGAP_LOG_INFO(
                "Calculating U_MM for multi-element component {}/{} ({} sparse points)",
                comp_idx + 1,
                new_multi_element_components.size(),
                comp->nSparsePoints()
            );
            auto U = U_MM(comp);
            const size_t m = U.nRows();
            assert(m == comp->nSparsePoints());

            for (size_t j = 0; j < m; j++) {
                for (size_t i = 0; i < m; i++) {
                    A(current_M + i, current_M + j) = U(i, j);
                }
            }
            current_M += m;
        }
    }
}
