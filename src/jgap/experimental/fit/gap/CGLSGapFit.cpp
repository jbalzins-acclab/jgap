#include "CGLSGapFit.hpp"

#include <atomic>
#include <cassert>
#include <limits>
#include <random>

#include "jgap/core/UnseqFor.hpp"
#include "jgap/core/linalg/Linalg.hpp"

namespace jgap {

    std::vector<Real> CGLSGapFit::findCoefficients(
        std::vector<ValuePtr<GapComponent>>& gap_components, const std::vector<Atoms>& training_data,
        std::vector<EnergyData>& energies_without_external, std::vector<Regularization>& sigmas_inverse
    ) {
        JGAP_LOG_INFO("Forming matrix A (RowMajor)");
        auto A = formMatrixA(gap_components, training_data, energies_without_external, sigmas_inverse);

        JGAP_LOG_INFO("Forming feature vector b");
        auto b = formVectorB(gap_components, energies_without_external, sigmas_inverse);

        JGAP_LOG_INFO("Doing linear algebra");
        auto c = leastSquares(A, b);

        return c;
    }

    std::vector<Real> CGLSGapFit::leastSquares(Matrix<RowMajor>& A, std::vector<Real>& b) {
        return linalg::solveLeastSquaresConjugateGradient(A, b, max_iterations, tolerance);
    }

    Matrix<RowMajor> CGLSGapFit::formMatrixA(
        const std::vector<ValuePtr<GapComponent>>& gap_components, const std::vector<Atoms>& training_data,
        const std::vector<EnergyData>& energy_data, const std::vector<Regularization>& sigmas_inverse
    ) const {
        struct StructureData {
            const Atoms& atoms;
            const EnergyData& energy_data;
            const Regularization& sigmas_inverse;
        };

        size_t r = 0;
        std::vector<std::pair<size_t, StructureData>> starting_rows_LK_NM;
        for (size_t i = 0; i < training_data.size(); i++) {
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
            "Forming in-memory {}x{}(~{}GB) A matrix (RowMajor)", r + c, c,
            static_cast<Real>((r + c) * c * sizeof(Real)) / 1024.0 / 1024.0 / 1024.0
        );
        Matrix<RowMajor> resulting_A(r + c, c);

        unseqForEach(
            starting_rows_LK_NM.begin(), starting_rows_LK_NM.end(),
            [&](const std::pair<size_t, StructureData>& starting_row_and_data) {
                const auto& [starting_row, data] = starting_row_and_data;

                fillInverseSigmaLK_NM(
                    gap_components, data.atoms, data.energy_data, data.sigmas_inverse, resulting_A, starting_row
                );
            }
        );

        unseqForEach(
            starting_points_K_MM.begin(), starting_points_K_MM.end(),
            [&](const std::array<size_t, 3>& rc_and_descriptor_id) {
                auto& [starting_row, starting_col, descriptor_id] = rc_and_descriptor_id;

                fillU_mm(starting_row, starting_col, gap_components[descriptor_id], resulting_A);
            }
        );

        return resulting_A;
    }

    std::vector<Real> CGLSGapFit::formVectorB(
        const std::vector<ValuePtr<GapComponent>>& components, const std::vector<EnergyData>& energy_data,
        const std::vector<Regularization>& sigmas_inverse
    ) {
        std::vector<Real> b;
        for (size_t i = 0; i < energy_data.size(); i++) {
            if (energy_data[i].energy.has_value()) {
                assert(sigmas_inverse[i].energy.has_value());
                b.push_back(energy_data[i].energy.value() * sigmas_inverse[i].energy.value());
            }

            if (energy_data[i].forces.has_value()) {
                assert(sigmas_inverse[i].forces.has_value());
                assert(sigmas_inverse[i].forces->size() == energy_data[i].forces->size());

                for (size_t j = 0; j < energy_data[i].forces->size(); j++) {
                    b.push_back(energy_data[i].forces->at(j).x * sigmas_inverse[i].forces->at(j).x);
                    b.push_back(energy_data[i].forces->at(j).y * sigmas_inverse[i].forces->at(j).y);
                    b.push_back(energy_data[i].forces->at(j).z * sigmas_inverse[i].forces->at(j).z);
                }
            }

            if (energy_data[i].virials.has_value()) {
                assert(sigmas_inverse[i].virials.has_value());
                b.push_back(energy_data[i].virials->xx * sigmas_inverse[i].virials->xx);
                b.push_back(energy_data[i].virials->xy * sigmas_inverse[i].virials->xy);
                b.push_back(energy_data[i].virials->xz * sigmas_inverse[i].virials->xz);
                b.push_back(energy_data[i].virials->yy * sigmas_inverse[i].virials->yy);
                b.push_back(energy_data[i].virials->yz * sigmas_inverse[i].virials->yz);
                b.push_back(energy_data[i].virials->zz * sigmas_inverse[i].virials->zz);
            }
        }

        for (const auto& component: components) {
            b.resize(b.size() + component->nSparsePoints(), 0.0_r);
        }

        return b;
    }

    void CGLSGapFit::fillInverseSigmaLK_NM(
        const std::vector<ValuePtr<GapComponent>>& gap_components, const Atoms& atoms, const EnergyData& energy_data,
        const Regularization& sigmas_inverse, Matrix<RowMajor>& A, const size_t starting_row
    ) {
        size_t contribution_column = 0;
        for (const auto& component: gap_components) {
            const auto max_cutoffs = component->getCutoffs();
            const NeighbourLists neighbour_list(atoms, max_cutoffs.maxOverall());

            const auto covariate_results = component->covariate(neighbour_list);

            if (!covariate_results.has_value()) {
                contribution_column += component->nSparsePoints();
                continue;
            }

            const auto& covariances = covariate_results.value();

            for (size_t sparse_idx = 0; sparse_idx < component->nSparsePoints(); sparse_idx++) {
                size_t currentRow = starting_row;

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
    }

    void CGLSGapFit::fillU_mm(
        const size_t starting_row, const size_t starting_col, const ValuePtr<GapComponent>& gap_component,
        Matrix<RowMajor>& A
    ) const {
        auto K_MM_block = gap_component->sparseToSparseCovariance();

        const size_t n = K_MM_block.nRows();
        for (size_t i = 0; i < n; i++) K_MM_block(i, i) += jitter;

        auto U_mm_block = choleskyDecomposition(K_MM_block);

        for (size_t i = 0; i < n; i++) {
            for (size_t j = 0; j < n; j++) {
                A(starting_row + i, starting_col + j) = U_mm_block(i, j);
            }
        }
    }

    Matrix<RowMajor> CGLSGapFit::choleskyDecomposition(Matrix<RowMajor>& matrix_block) {
        return linalg::choleskyDecomposition<MatrixLayout::RowMajor>(matrix_block);
    }
} // namespace jgap
