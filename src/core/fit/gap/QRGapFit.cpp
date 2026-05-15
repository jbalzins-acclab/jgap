#include "QRGapFit.hpp"

#include <random>

#include "io/log/StdoutLogger.hpp"
#include "utils/Utils.hpp"

#include <tbb/parallel_for_each.h>
#include <cassert>

#include <Eigen/Dense>

namespace jgap {

    QRGapFit::QRGapFit(const Real jitter) : jitter(jitter) {
    }

    std::vector<Real> QRGapFit::mainFit(std::vector<IGapComponent::Ptr> &gap_components,
                                        std::vector<NeighbourList> &neighbour_lists,
                                        std::vector<EnergyData> energies_without_external,
                                        std::vector<Regularization> &sigmas_inverse) {

        JGAP_LOG_INFO("Forming matrix A");
        auto A = formMatrixA(gap_components, neighbour_lists, energies_without_external, sigmas_inverse);

        JGAP_LOG_INFO("Forming feature vector b");
        auto b = formVectorB(gap_components, energies_without_external, sigmas_inverse);

        JGAP_LOG_INFO("Doing linear algebra");
        auto c = leastSquares(A, b);

        return c;
    }

    std::vector<Real> QRGapFit::leastSquares(EigenMatrix &A, EigenVector &b) {
        JGAP_LOG_DEBUG("Init Eigen::HouseholderQR");
        const Eigen::HouseholderQR<Eigen::Ref<EigenMatrix>> qr(A);

        JGAP_LOG_DEBUG("Q^t");
        auto Qt = qr.householderQ().transpose();

        JGAP_LOG_DEBUG("Q^t * b");
        auto Qt_b =  Qt * b;

        JGAP_LOG_DEBUG("R");
        auto R = qr.matrixQR().topLeftCorner(A.cols(), A.cols());

        JGAP_LOG_DEBUG("R^-1 * Q^t_b");
        EigenVector c = R.triangularView<Eigen::Upper>().solve(Qt_b.head(A.cols()));

        return std::vector<Real>{c.data(), c.data() + c.size()};
    }

    EigenMatrix QRGapFit::formMatrixA(const std::vector<IGapComponent::Ptr> &gap_components,
                                      const std::vector<NeighbourList> &neighbour_lists,
                                      const std::vector<EnergyData> &energy_data,
                                      const std::vector<Regularization> &sigmas_inverse) const {

        struct StructureData {
            const NeighbourList& neighbour_lists;
            const EnergyData& energy_data;
            const Regularization& sigmas_inverse;
        };

        size_t r = 0;
        std::vector<std::pair<size_t, StructureData>> startingRowsK_nm;
        for (int i = 0; i < neighbour_lists.size(); i++) {
            startingRowsK_nm.emplace_back(r, StructureData{
                neighbour_lists[i],
                energy_data[i],
                sigmas_inverse[i]
            });
            if (energy_data[i].energy.has_value()) r += 1;
            if (energy_data[i].forces.has_value()) r += 3 * energy_data[i].forces->size();
            if (energy_data[i].virials.has_value()) r += 6;
        }

        std::vector<std::array<size_t, 3>> startingPointsK_mm;
        size_t c = 0;
        for (size_t i = 0; i < gap_components.size(); i++) {
            startingPointsK_mm.push_back({r + c, c, i});
            c += gap_components[i]->nSparsePoints();
        }

        JGAP_LOG_INFO("Forming in-memory {}x{}(~{}GB) A matrix",
                      r+c, c, (r+c)*c * sizeof(double) / 1024.0 / 1024.0 / 1024.0);
        EigenMatrix resultingA = EigenMatrix::Zero(r + c, c);

        std::atomic counter(0);
        tbb::parallel_for_each(
            startingRowsK_nm.begin(), startingRowsK_nm.end(),
            [&](const std::pair<size_t, StructureData>& structId) {

                auto &[starting_row, struct_data] = structId;

                size_t progress = ++counter;
                if (progress % std::max(startingRowsK_nm.size() / 100, 1uz) == 0) {
                    JGAP_LOG_INFO(
                        "K_nm matrix formation progress: {} of {} ({}%)",
                        progress, startingRowsK_nm.size(), progress * 100 / startingRowsK_nm.size()
                    );
                }

                fillInverseSigmaK_nm(
                    gap_components,
                    struct_data.neighbour_lists,
                    struct_data.energy_data,
                    struct_data.sigmas_inverse,
                    resultingA,
                    starting_row
                    );
            }
        );

        tbb::parallel_for_each(
            startingPointsK_mm.begin(), startingPointsK_mm.end(),
            [&](const std::array<size_t, 3>& rcAndDescriptorId) {
                auto &[starting_row, starting_col, descriptorId] = rcAndDescriptorId;

                JGAP_LOG_INFO("K_mm for descriptor {}", descriptorId);
                fillU_mm(
                    starting_row,
                    starting_col,
                    gap_components[descriptorId],
                    resultingA
                    );
            }
        );

        return resultingA;
    }

    EigenVector QRGapFit::formVectorB(const std::vector<IGapComponent::Ptr> &components,
                                      const std::vector<EnergyData> &energy_data,
                                      const std::vector<Regularization> &sigmas_inverse) {

        std::vector<Real> b;
        for (size_t i = 0; i < energy_data.size(); i++) {
            if (energy_data[i].energy.has_value()) {
                assert(sigmas_inverse[i].energy.has_value());
                b.push_back(energy_data[i].energy.value() * sigmas_inverse[i].energy.value());
            }

            if (energy_data[i].forces.has_value()) {
                assert(sigmas_inverse[i].forces.has_value());
                assert(sigmas_inverse[i].forces->size() == energy_data[i].forces->size());

                for (int j = 0; j < energy_data[i].forces->size(); j++) {
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
            b.resize(b.size() + component->nSparsePoints(), 0.0);
        }

        return Eigen::Map<EigenVector>(b.data(), b.size());
    }

    void QRGapFit::fillInverseSigmaK_nm(const std::vector<IGapComponent::Ptr> &gap_components,
                                        const NeighbourList& neighbour_list,
                                        const EnergyData& energy_data,
                                        const Regularization& sigmas_inverse,
                                        EigenMatrix &A,
                                        size_t starting_row) {

        size_t contributionColumn = 0;
        for (const auto& gap_component: gap_components) {
            auto contributions = gap_component->covariate(neighbour_list);

            for (const auto& contribution: contributions) {

                size_t currentRow = starting_row;

                if (energy_data.energy.has_value()) {
                    A(currentRow++, contributionColumn) = contribution.value * sigmas_inverse.energy.value();
                }

                if (energy_data.forces.has_value() && !contribution.forces.empty()) {
                    for (size_t rowInc = 0; rowInc < contribution.forces.size(); rowInc++) {

                        const auto force = contribution.forces[rowInc];
                        const auto fSigmasInverse = (*sigmas_inverse.forces)[rowInc];

                        A(currentRow++, contributionColumn) = force.x * fSigmasInverse.x;
                        A(currentRow++, contributionColumn) = force.y * fSigmasInverse.y;
                        A(currentRow++, contributionColumn) = force.z * fSigmasInverse.z;
                    }
                }

                if (energy_data.virials.has_value()) {
                    auto [xx, xy, xz, yy, yz, zz] = contribution.virials;
                    auto [sigma_xx, sigma_xy, sigma_xz, sigma_yy, sigma_yz, sigma_zz]
                        = sigmas_inverse.virials.value();

                    A(currentRow++, contributionColumn) = xx * sigma_xx;
                    A(currentRow++, contributionColumn) = xy * sigma_xy;
                    A(currentRow++, contributionColumn) = xz * sigma_xz;

                    A(currentRow++, contributionColumn) = yy * sigma_yy;
                    A(currentRow++, contributionColumn) = yz * sigma_yz;

                    A(currentRow++, contributionColumn) = zz * sigma_zz;
                }

                contributionColumn++;
            }
        }
    }

    void QRGapFit::fillU_mm(const size_t startingRow,
                            const size_t startingCol,
                            const IGapComponent::Ptr& gap_component,
                            EigenMatrix &A) const {

        size_t nPreviousKernels = 0;
        for (auto &K_mmBlock: gap_component->sparseToSparseCovariance()) {

            const size_t n = K_mmBlock.rows();
            for (size_t i = 0; i < n; i++) K_mmBlock(i, i) += jitter;

            auto K_mmConverted = convertToEigen(K_mmBlock);
            auto U_mmBlock = choleskyDecomposition(K_mmConverted);

            for (size_t i = 0; i < n; i++) {
                for (size_t j = 0; j < n; j++) {
                    A(startingRow + nPreviousKernels + i, startingCol + nPreviousKernels + j) = U_mmBlock(i, j);
                }
            }

            nPreviousKernels += n;
        }
    }

    EigenMatrix QRGapFit::choleskyDecomposition(EigenMatrix &matrix_block) {
        Eigen::LLT<EigenMatrix> llt(matrix_block);
        if (llt.info() != Eigen::Success) {
            JGAP_LOG_ERROR("Cholesky decomposition failed: matrix not positive definite", true);
        }

        return llt.matrixU();
    }

    EigenMatrix QRGapFit::convertToEigen(Matrix &matrix_block) {
        return Eigen::Map<Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>(
            matrix_block.rawData().data(), matrix_block.rows(), matrix_block.columns()
            );
    }
}
