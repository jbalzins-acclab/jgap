#include "IterativeGapFit.hpp"

#include <Eigen/Dense>
#include <Eigen/IterativeLinearSolvers>
#include <atomic>
#include <cassert>
#include <limits>
#include <random>

#include "jgap/core/UnseqFor.hpp"

namespace jgap {

    using EigenMatrixCol = Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>;
    using EigenMatrixRow = Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
    using EigenVectorCol = Eigen::Matrix<Real, Eigen::Dynamic, 1>;

    std::vector<Real> IterativeGapFit::findCoefficients(
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

    std::vector<Real> IterativeGapFit::leastSquares(Matrix<RowMajor>& A, std::vector<Real>& b) {
        const size_t R = A.nRows();
        const size_t C = A.nColumns();

        Eigen::Map<EigenMatrixRow> A_map(A.flatData().data(), R, C);
        Eigen::Map<EigenVectorCol> b_map(b.data(), b.size());

        JGAP_LOG_INFO("Solving least squares via LeastSquaresConjugateGradient directly on A ({}x{})", R, C);

        Eigen::LeastSquaresConjugateGradient<EigenMatrixRow> lscg;
        if (max_iterations.has_value()) {
            lscg.setMaxIterations(max_iterations.value());
        }
        if (tolerance.has_value()) {
            lscg.setTolerance(tolerance.value());
        }

        lscg.compute(A_map);
        if (lscg.info() != Eigen::Success) {
            JGAP_LOG_ERROR("LeastSquaresConjugateGradient setup failed", true);
        }

        EigenVectorCol c = lscg.solve(b_map);

        JGAP_LOG_INFO("LSCG finished: iterations = {}, relative error = {}", lscg.iterations(), lscg.error());

        return std::vector<Real>{c.data(), c.data() + c.size()};
    }

    Matrix<RowMajor> IterativeGapFit::formMatrixA(
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
            (r + c) * c * sizeof(double) / 1024.0 / 1024.0 / 1024.0
        );
        Matrix<RowMajor> resulting_A(r + c, c);

        std::atomic<size_t> counter(0);
        unseqForEach(
            starting_rows_LK_NM.begin(), starting_rows_LK_NM.end(),
            [&](const std::pair<size_t, StructureData>& structId) {
                auto& [starting_row, struct_data] = structId;

                size_t progress = ++counter;
                if (progress % std::max(starting_rows_LK_NM.size() / 100, 1uz) == 0) {
                    JGAP_LOG_INFO(
                        "LK_NM matrix formation progress: {} of {} ({}%)", progress, starting_rows_LK_NM.size(),
                        progress * 100 / starting_rows_LK_NM.size()
                    );
                }

                fillInverseSigmaLK_NM(
                    gap_components, struct_data.atoms, struct_data.energy_data, struct_data.sigmas_inverse, resulting_A,
                    starting_row
                );
            }
        );

        unseqForEach(
            starting_points_K_MM.begin(), starting_points_K_MM.end(),
            [&](const std::array<size_t, 3>& rc_and_descriptor_id) {
                auto& [starting_row, starting_col, descriptor_id] = rc_and_descriptor_id;

                JGAP_LOG_INFO("K_MM for descriptor {}", descriptor_id);
                fillU_mm(starting_row, starting_col, gap_components[descriptor_id], resulting_A);
            }
        );

        return resulting_A;
    }

    std::vector<Real> IterativeGapFit::formVectorB(
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

        for (auto& component: components) {
            b.resize(b.size() + component->nSparsePoints(), 0.0_r);
        }

        return b;
    }

    void IterativeGapFit::fillInverseSigmaLK_NM(
        const std::vector<ValuePtr<GapComponent>>& gap_components, const Atoms& atoms, const EnergyData& energy_data,
        const Regularization& sigmas_inverse, Matrix<RowMajor>& A, size_t starting_row
    ) {
        std::map<Real, NeighbourLists> neighbour_lists;
        for (const auto& gap_component: gap_components) {
            Real cutoff = gap_component->getCutoff();
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

    void IterativeGapFit::fillU_mm(
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

    Matrix<RowMajor> IterativeGapFit::choleskyDecomposition(Matrix<RowMajor>& matrix_block) {
        using EigenMatrixRow = Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
        Eigen::Map<EigenMatrixRow> map(matrix_block.flatData().data(), matrix_block.nRows(), matrix_block.nColumns());

        Eigen::LLT<Eigen::Map<EigenMatrixRow>> llt(map);
        if (llt.info() != Eigen::Success) {
            JGAP_LOG_ERROR("Cholesky decomposition failed: matrix not positive definite", true);
        }

        Matrix<RowMajor> result(matrix_block.nRows(), matrix_block.nColumns());
        Eigen::Map<EigenMatrixRow> res_map(result.flatData().data(), result.nRows(), result.nColumns());
        res_map = llt.matrixU();
        return result;
    }
}
