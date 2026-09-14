#include "QRGapFit.hpp"

#include <atomic>
#include <cassert>
#include <cmath>

#include "jgap/core/UnseqFor.hpp"
#include "jgap/core/linalg/Linalg.hpp"

namespace jgap {

    void QRGapFit::findCoefficients(
        GapPotential& to_be_fit,
        const std::vector<Atoms>& training_data,
        std::vector<EnergyData>& energies_without_external,
        std::vector<Regularization>& sigmas_inverse
    ) {
        JGAP_LOG_INFO("Forming matrix A");
        auto A = formAugmentedCovarianceMatrixA(
            to_be_fit.components,
            training_data,
            energies_without_external,
            sigmas_inverse
            );

        JGAP_LOG_INFO("Forming feature vector b");
        auto b = formNormalizedAugmentedTargetVectorB(to_be_fit.components, energies_without_external, sigmas_inverse);

        JGAP_LOG_INFO("Doing linear algebra");
        auto c = leastSquares(A, b);

        to_be_fit.setCoefficients(c);
    }

    std::vector<double> QRGapFit::leastSquares(Matrix<ColumnMajor>& A, std::vector<double>& b) {
        return linalg::solveLeastSquaresHouseholderQR(A, b);
    }

    Matrix<ColumnMajor> QRGapFit::formAugmentedCovarianceMatrixA(
        const std::vector<ValuePtr<GapComponent>>& gap_components,
        const std::vector<Atoms>& training_data,
        const std::vector<EnergyData>& energy_data,
        const std::vector<Regularization>& sigmas_inverse
    ) const {
        struct StructureData {
            const Atoms& atoms;
            const EnergyData& energy_data;
            const Regularization& sigmas_inverse;
        };

        size_t r = 0;
        std::vector<std::pair<size_t, StructureData>> starting_rows_LK_NM;
        for (int i = 0; i < training_data.size(); i++) {
            starting_rows_LK_NM.emplace_back(r, StructureData{training_data[i], energy_data[i], sigmas_inverse[i]});
            if (energy_data[i].energy.has_value()) r += 1;
            if (energy_data[i].forces.has_value()) r += 3 * energy_data[i].forces->size();
            if (energy_data[i].virials.has_value()) r += 6;
        }

        std::vector<std::array<size_t, 3>> starting_points_K_MM;
        size_t c = 0;
        for (size_t i = 0; i < gap_components.size(); i++) {
            starting_points_K_MM.push_back({r + c, c, i});
            c += gap_components[i]->nSparsePoints();
        }

        JGAP_LOG_INFO(
            "Forming in-memory {}x{}(~{}GB) A matrix", r + c, c, (r + c) * c * sizeof(double) / 1024.0 / 1024.0 / 1024.0
        );
        Matrix<ColumnMajor> resulting_A(r + c, c);

        std::atomic counter(0);
        unseqForEach(
            starting_rows_LK_NM.begin(),
            starting_rows_LK_NM.end(),
            [&](const std::pair<size_t, StructureData>& structId) {
                auto& [starting_row, struct_data] = structId;

                size_t progress = ++counter;
                if (progress % std::max(starting_rows_LK_NM.size() / 100, 1uz) == 0) {
                    JGAP_LOG_INFO(
                        "LK_NM matrix formation progress: {} of {} ({}%)",
                        progress,
                        starting_rows_LK_NM.size(),
                        progress * 100 / starting_rows_LK_NM.size()
                    );
                }

                auto A_entry = formInverseSigmaLK_NMForEntry(
                    gap_components,
                    struct_data.atoms,
                    struct_data.energy_data,
                    struct_data.sigmas_inverse
                );
                for (size_t col = 0; col < c; ++col) {
                    for (size_t row = 0; row < A_entry.nRows(); ++row) {
                        resulting_A(starting_row + row, col) = A_entry(row, col);
                    }
                }
            }
        );

        unseqForEach(
            starting_points_K_MM.begin(),
            starting_points_K_MM.end(),
            [&](const std::array<size_t, 3>& rc_and_descriptor_id) {
                auto& [starting_row, starting_col, descriptor_id] = rc_and_descriptor_id;

                JGAP_LOG_INFO("U_MM for descriptor {}", descriptor_id);
                auto U = U_MM(gap_components[descriptor_id]);
                const size_t n = U.nRows();
                for (size_t i = 0; i < n; i++) {
                    for (size_t j = 0; j < n; j++) {
                        resulting_A(starting_row + i, starting_col + j) = U(i, j);
                    }
                }
            }
        );

        return resulting_A;
    }

    std::vector<double> QRGapFit::formTargetVectorBForEntry(
        const EnergyData& energy_data,
        const Regularization& sigmas_inverse
    ) {
        std::vector<double> b;
        if (energy_data.energy.has_value()) {
            assert(sigmas_inverse.energy.has_value());
            b.push_back(energy_data.energy.value() * sigmas_inverse.energy.value());
        }

        if (energy_data.forces.has_value()) {
            assert(sigmas_inverse.forces.has_value());
            assert(sigmas_inverse.forces->size() == energy_data.forces->size());

            for (size_t j = 0; j < energy_data.forces->size(); j++) {
                b.push_back(energy_data.forces->at(j).x * sigmas_inverse.forces->at(j).x);
                b.push_back(energy_data.forces->at(j).y * sigmas_inverse.forces->at(j).y);
                b.push_back(energy_data.forces->at(j).z * sigmas_inverse.forces->at(j).z);
            }
        }

        if (energy_data.virials.has_value()) {
            assert(sigmas_inverse.virials.has_value());
            b.push_back(energy_data.virials->xx * sigmas_inverse.virials->xx);
            b.push_back(energy_data.virials->xy * sigmas_inverse.virials->xy);
            b.push_back(energy_data.virials->xz * sigmas_inverse.virials->xz);
            b.push_back(energy_data.virials->yy * sigmas_inverse.virials->yy);
            b.push_back(energy_data.virials->yz * sigmas_inverse.virials->yz);
            b.push_back(energy_data.virials->zz * sigmas_inverse.virials->zz);
        }

        return b;
    }

    std::vector<double> QRGapFit::formNormalizedAugmentedTargetVectorB(
        const std::vector<ValuePtr<GapComponent>>& components,
        const std::vector<EnergyData>& energy_data,
        const std::vector<Regularization>& sigmas_inverse
    ) {
        std::vector<double> b;
        for (size_t i = 0; i < energy_data.size(); i++) {
            auto b_chunk = formTargetVectorBForEntry(energy_data[i], sigmas_inverse[i]);
            b.insert(b.end(), b_chunk.begin(), b_chunk.end());
        }

        for (auto& component: components) {
            b.resize(b.size() + component->nSparsePoints(), 0.0);
        }

        return b;
    }

    Matrix<ColumnMajor> QRGapFit::formInverseSigmaLK_NMForEntry(
        const std::vector<ValuePtr<GapComponent>>& gap_components,
        const Atoms& atoms,
        const EnergyData& energy_data,
        const Regularization& sigmas_inverse
    ) {
        size_t n_rows = 0;
        if (energy_data.energy.has_value()) n_rows += 1;
        if (energy_data.forces.has_value()) n_rows += 3 * atoms.nAtoms();
        if (energy_data.virials.has_value()) n_rows += 6;

        size_t n_cols = 0;
        for (const auto& comp: gap_components) {
            n_cols += comp->nSparsePoints();
        }

        Matrix<ColumnMajor> A(n_rows, n_cols);

        std::map<double, NeighbourLists> neighbour_lists;
        for (const auto& gap_component: gap_components) {
            double cutoff = gap_component->getCutoff();
            if (!neighbour_lists.contains(cutoff)) {
                neighbour_lists.insert({cutoff, NeighbourLists(atoms, cutoff)});
            }
        }

        size_t contribution_column = 0;
        for (const auto& gap_component: gap_components) {
            auto& neighbour_list = neighbour_lists.at(gap_component->getCutoff());
            auto covariances_opt = gap_component->covariate(neighbour_list);

            if (!covariances_opt.has_value()) {
                contribution_column += gap_component->nSparsePoints();
                continue;
            }

            auto& covariances = covariances_opt.value();

            for (size_t sparse_idx = 0; sparse_idx < gap_component->nSparsePoints(); sparse_idx++) {
                size_t currentRow = 0;

                if (energy_data.energy.has_value()) {
                    A(currentRow++, contribution_column) =
                        covariances.energy(sparse_idx) * sigmas_inverse.energy.value();
                }

                if (energy_data.forces.has_value()) {
                    for (size_t rowInc = 0; rowInc < atoms.nAtoms(); rowInc++) {
                        const auto& force = covariances.force(sparse_idx, rowInc);
                        const auto fSigmasInverse = (*sigmas_inverse.forces)[rowInc];

                        A(currentRow++, contribution_column) = force.x * fSigmasInverse.x;
                        A(currentRow++, contribution_column) = force.y * fSigmasInverse.y;
                        A(currentRow++, contribution_column) = force.z * fSigmasInverse.z;
                    }
                }

                if (energy_data.virials.has_value()) {
                    auto [xx, xy, xz, yy, yz, zz] = covariances.virials(sparse_idx);
                    auto [sigma_xx, sigma_xy, sigma_xz, sigma_yy, sigma_yz, sigma_zz] = sigmas_inverse.virials.value();

                    A(currentRow++, contribution_column) = xx * sigma_xx;
                    A(currentRow++, contribution_column) = xy * sigma_xy;
                    A(currentRow++, contribution_column) = xz * sigma_xz;

                    A(currentRow++, contribution_column) = yy * sigma_yy;
                    A(currentRow++, contribution_column) = yz * sigma_yz;

                    A(currentRow++, contribution_column) = zz * sigma_zz;
                }

                contribution_column++;
            }
        }

        return A;
    }

    Matrix<ColumnMajor> QRGapFit::U_MM(const ValuePtr<GapComponent>& gap_component) const {
        auto K_MM_block = gap_component->K_MM();

        const size_t n = K_MM_block.nRows();
        for (size_t i = 0; i < n; i++) K_MM_block(i, i) += jitter;

        return choleskyDecomposition(K_MM_block);
    }

    Matrix<ColumnMajor> QRGapFit::choleskyDecomposition(Matrix<RowMajor>& matrix_block) {
        return linalg::choleskyDecomposition<MatrixLayout::ColumnMajor>(matrix_block);
    }
}
