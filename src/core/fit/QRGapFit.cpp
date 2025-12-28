#include "core/fit/QRGapFit.hpp"

#include <random>

#include "core/neighbours/NeighbourFinder.hpp"
#include "io/log/StdoutLogger.hpp"
#include "utils/Utils.hpp"
#include "core/matrices/regularization/RegularizationRules.hpp"

#include <tbb/parallel_for_each.h>

#include <Eigen/Sparse>

namespace jgap {

    std::shared_ptr<QRGapFit> QRGapFit::fromDataNode(const DataNode& params) {

        std::map<std::string, std::shared_ptr<Descriptor>> descriptors{};
        for (const auto& [label, descriptorParams]: require(params, "descriptors").asObject()) {
            descriptors[label] = REGISTRY_GET(Descriptor, descriptorParams);
        }

        auto regularizationRules = REGISTRY_GET(RegularizationRules, require(params, "regularization_rules"));
        double jitter = params.getOrDefault("jitter", 1e-8);

        return std::make_shared<QRGapFit>(descriptors, regularizationRules, jitter);
    }

    QRGapFit::QRGapFit(
        const std::map<std::string, std::shared_ptr<Descriptor>> &descriptors,
        const std::shared_ptr<RegularizationRules> &regularization_rules,
        const double jitter
        ) : descriptors_(descriptors), regularization_rules_(regularization_rules), _jitter(jitter) {
    }

    std::shared_ptr<Potential> QRGapFit::fit(const std::vector<AtomicStructure>& training_data) {
        JGAP_LOG_INFO("Starting JGAP fit");

        std::vector<AtomicStructure> _trainingData(training_data);
        JGAP_LOG_INFO("Checking sparse points");

        double maxCutoff = 0;
        auto descriptorsAsVec = std::vector<std::shared_ptr<Descriptor>>();
        for (const auto& descriptor: descriptors_ | std::views::values) {
            descriptorsAsVec.push_back(descriptor);

            NeighbourFinder::findNeighbours(_trainingData, descriptor->getCutoff().maxOverall());
            descriptor->setupSparseKernels(_trainingData);

            maxCutoff = std::max(maxCutoff, descriptor->getCutoff().maxOverall());
        }

        JGAP_LOG_DEBUG("Full neighbour-std::list");
        NeighbourFinder::findNeighbours(_trainingData, maxCutoff);

        //// ----------------------------------------------------------------------------------------------------

        JGAP_LOG_INFO("Setting up regularization parameters");
        for (auto& structure: _trainingData) {
            regularization_rules_->fillSigmas(structure);
        }

        //// ----------------------------------------------------------------------------------------------------

        JGAP_LOG_INFO("Making matrix A");
        auto A = makeA(descriptorsAsVec, _trainingData);

        JGAP_LOG_INFO("Making feature std::vector b");
        auto b = makeB(descriptorsAsVec, _trainingData);

        //// ----------------------------------------------------------------------------------------------------

        JGAP_LOG_INFO("Doing linear algebra");
        auto c = leastSquares(A, b);

        size_t counter = 0;
        for (const auto &descriptor: descriptors_ | std::views::values) {
            for (const auto& kernel: descriptor->getKernels()) {
                kernel->coefficient = c[counter++];
            }
        }

        return std::make_shared<GapPotential>(descriptors_);
    }

    std::vector<double> QRGapFit::leastSquares(Eigen::MatrixXd &A, Eigen::VectorXd &b) {
        JGAP_LOG_DEBUG("Init Eigen::HouseholderQR");
        const Eigen::HouseholderQR<Eigen::Ref<Eigen::MatrixXd>> qr(A);

        JGAP_LOG_DEBUG("Q^t");
        auto Qt = qr.householderQ().transpose();

        JGAP_LOG_DEBUG("Q^t * b");
        auto Qt_b =  Qt * b;

        JGAP_LOG_DEBUG("R");
        auto R = qr.matrixQR().topLeftCorner(A.cols(), A.cols());

        JGAP_LOG_DEBUG("R^-1 * Q^t_b");
        Eigen::VectorXd c = R.triangularView<Eigen::Upper>().solve(Qt_b.head(A.cols()));

        return std::vector<double>{c.data(), c.data() + c.size()};
    }

    Eigen::MatrixXd QRGapFit::makeA(const std::vector<std::shared_ptr<Descriptor>> &descriptors,
                                    const std::vector<AtomicStructure> &atomic_structures) const {

        size_t r = 0;
        std::vector<std::pair<size_t, AtomicStructure>> startingRowsK_nm;
        for (const auto& structure : atomic_structures) {
            startingRowsK_nm.emplace_back(r, structure);
            if (structure.energy.has_value()) r += 1;
            if (structure.forces.has_value()) r += 3 * structure.size();
            if (structure.virials.has_value()) r += 6;
        }

        std::vector<std::array<size_t, 3>> startingPointsK_mm;
        size_t c = 0;
        for (size_t i = 0; i < descriptors.size(); i++) {
            startingPointsK_mm.push_back({r + c, c, i});
            auto kernels = descriptors[i]->getKernels();
            c += descriptors[i]->getKernels().size();
            JGAP_LOG_INFO("nkernels = {}", descriptors[i]->getKernels().size());
        }

        JGAP_LOG_DEBUG("Forming in-memory {}x{}(~{}GB) A matrix",
                      r+c, c, (r+c)*c * sizeof(double) / 1024.0 / 1024.0 / 1024.0);
        Eigen::MatrixXd resultingA = Eigen::MatrixXd::Zero(r + c, c);

        std::atomic counter(0);
        tbb::parallel_for_each(
            startingRowsK_nm.begin(), startingRowsK_nm.end(),
            [&](const std::pair<size_t, AtomicStructure>& structId) {

                size_t progress = ++counter;
                if (progress % std::max(startingRowsK_nm.size() / 100, 1uz) == 0) {
                    JGAP_LOG_DEBUG(
                        "K_nm matrix formation progress: {} of {} ({}%)",
                        progress, startingRowsK_nm.size(), progress * 100 / startingRowsK_nm.size()
                    );
                }

                fillInverseSigmaK_nm(descriptors, structId.second, resultingA, structId.first);
            }
        );

        tbb::parallel_for_each(
            startingPointsK_mm.begin(), startingPointsK_mm.end(),
            [&](const std::array<size_t, 3>& descriptorId) {
                JGAP_LOG_DEBUG("K_mm for descriptor {}", descriptorId[2]);
                fillU_mm(descriptorId[0], descriptorId[1], *descriptors[descriptorId[2]], resultingA);
            }
        );

        return resultingA;
    }

    Eigen::VectorXd QRGapFit::makeB(const std::vector<std::shared_ptr<Descriptor>> &descriptors,
                                    const std::vector<AtomicStructure> &atomic_structures) {
        std::vector<double> b;
        for (auto& structure: atomic_structures) {
            if (structure.energy.has_value()) {
                b.push_back(structure.energy.value() * structure.energySigmaInverse.value());
            }
            if (structure.forces.has_value()) {
                for (const auto& atom: structure) {
                    b.push_back(atom.force().x * atom.forceSigmasInverse().x);
                    b.push_back(atom.force().y * atom.forceSigmasInverse().y);
                    b.push_back(atom.force().z * atom.forceSigmasInverse().z);
                }
            }
            if (structure.virials.has_value()) {
                b.push_back(structure.virials.value()[0].x * structure.virialSigmasInverse.value()[0].x);
                b.push_back(structure.virials.value()[0].y * structure.virialSigmasInverse.value()[0].y);
                b.push_back(structure.virials.value()[0].z * structure.virialSigmasInverse.value()[0].z);

                b.push_back(structure.virials.value()[1].y * structure.virialSigmasInverse.value()[1].y);
                b.push_back(structure.virials.value()[1].z * structure.virialSigmasInverse.value()[1].z);

                b.push_back(structure.virials.value()[2].z * structure.virialSigmasInverse.value()[2].z);
            }
        }

        for (auto& descriptor: descriptors) {
            b.resize(b.size() + descriptor->getKernels().size());
        }

        return Eigen::Map<Eigen::VectorXd>(b.data(), b.size());
    }

    void QRGapFit::fillInverseSigmaK_nm(const std::vector<std::shared_ptr<Descriptor>> &descriptors,
                                        const AtomicStructure &atomic_structure,
                                        Eigen::MatrixXd &A,
                                        const size_t starting_row) {
        size_t contributionColumn = 0;
        for (const auto& descriptor : descriptors) {
            auto contributions = descriptor->covariate(atomic_structure);

            for (const auto& contribution: contributions) {

                size_t currentRow = starting_row;

                if (atomic_structure.energy.has_value()) {
                    A(currentRow++, contributionColumn) =
                        contribution.total * atomic_structure.energySigmaInverse.value();
                }

                if (atomic_structure.forces.has_value() && !contribution.forces.empty()) {
                    for (size_t rowInc = 0; rowInc < contribution.forces.size(); rowInc++) {

                        const auto force = contribution.forces[rowInc]; // FORCE IS NEGATIVE
                        const auto fSigmasInverse = (*atomic_structure.forceSigmasInverse)[rowInc];

                        A(currentRow++, contributionColumn) = force.x * fSigmasInverse.x;
                        A(currentRow++, contributionColumn) = force.y * fSigmasInverse.y;
                        A(currentRow++, contributionColumn) = force.z * fSigmasInverse.z;
                    }
                }

                if (atomic_structure.virials.has_value()) {
                    auto virials = contribution.virials;
                    A(currentRow++, contributionColumn) =
                        virials[0].x * atomic_structure.virialSigmasInverse.value()[0].x;
                    A(currentRow++, contributionColumn) =
                        virials[0].y * atomic_structure.virialSigmasInverse.value()[0].y;
                    A(currentRow++, contributionColumn) =
                        virials[0].z * atomic_structure.virialSigmasInverse.value()[0].z;

                    A(currentRow++, contributionColumn) =
                        virials[1].y * atomic_structure.virialSigmasInverse.value()[1].y;
                    A(currentRow++, contributionColumn) =
                        virials[1].z * atomic_structure.virialSigmasInverse.value()[1].z;

                    A(currentRow++, contributionColumn) =
                        virials[2].z * atomic_structure.virialSigmasInverse.value()[2].z;
                }

                contributionColumn++;
            }
        }
    }

    void QRGapFit::fillU_mm(const size_t startingRow, const size_t startingCol,
                                Descriptor &descriptor, Eigen::MatrixXd &A) const {

        size_t nPreviousKernels = 0;
        for (auto &K_mmBlock: descriptor.selfCovariate()) {

            const size_t n = K_mmBlock->rows();
            for (size_t i = 0; i < n; i++) (*K_mmBlock)(i, i) += _jitter;

            auto K_mmConverted = convertToEigen(*K_mmBlock);
            auto U_mmBlock = choleskyDecomposition(K_mmConverted);

            for (size_t i = 0; i < n; i++) {
                for (size_t j = 0; j < n; j++) {
                    A(startingRow + nPreviousKernels + i, startingCol + nPreviousKernels + j) = U_mmBlock(i, j);
                }
            }

            nPreviousKernels += n;
        }
    }

    Eigen::MatrixXd QRGapFit::choleskyDecomposition(Eigen::MatrixXd &matrix_block) {
        Eigen::LLT<Eigen::MatrixXd> llt(matrix_block);
        if (llt.info() != Eigen::Success) {
            JGAP_LOG_ERROR("Cholesky decomposition failed: matrix not positive definite", true);
        }

        return llt.matrixU();
    }

    Eigen::MatrixXd QRGapFit::convertToEigen(MatrixBlock &matrix_block) {
        return Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>(
            matrix_block.rawData().data(), matrix_block.rows(), matrix_block.columns()
            );
    }
}
