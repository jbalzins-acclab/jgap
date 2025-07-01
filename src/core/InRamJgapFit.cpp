#include "core/InRamJgapFit.hpp"

#include "core/matrices/sigmas/SimpleSigmaRules.hpp"
#include "core/neighbours/NeighbourFinder.hpp"
#include "io/log/StdoutLogger.hpp"
#include "utils/Utils.hpp"

#include <tbb/parallel_for_each.h>

namespace jgap {

    InRamJgapFit::InRamJgapFit(const vector<shared_ptr<Descriptor>>& descriptors,
                     const vector<shared_ptr<Potential>>& externalPotentials,
                     const shared_ptr<SigmaRules>& sigmaRules,
                     const double jitter) : _jitter(jitter) {

        StdoutLogger::initIfNotInitialized();

        _externalPotentials = externalPotentials;
        _descriptors = descriptors;

        if (sigmaRules != nullptr) {
            _sigmaRules = sigmaRules;
        } else {
            _sigmaRules = make_shared<SimpleSigmaRules>(0.001, 0.04, 5, 5);
        }
    }

    shared_ptr<JgapPotential> InRamJgapFit::fitGAP(const vector<AtomicStructure>& trainingData) const {
        Logger::logger->info("Starting JGAP fit");

        const auto maxCutoffDescriptor = *ranges::max_element(_descriptors.begin(), _descriptors.end(),
                [](const shared_ptr<Descriptor> &a, const shared_ptr<Descriptor> &b) {
                    return a->getCutoff() < b->getCutoff();
                });
        auto dataToBeFit = trainingData;
        NeighbourFinder::findNeighbours(dataToBeFit, maxCutoffDescriptor -> getCutoff());
        dataToBeFit = subtractExternalContributions(dataToBeFit);

        Logger::logger->info("Checking sparse points");
        for (const auto& descriptor: _descriptors) {
            if (descriptor->nSparsePoints() == 0) {
                Logger::logger->info(format("No sparse points found in descriptor {} -> setting from data",
                                        descriptor->serialize().dump()));
                descriptor->setSparsePoints(dataToBeFit);
            }
        }

        for (auto& structure: dataToBeFit) {
            _sigmaRules->fillSigmas(structure);
        }

        //// ----------------------------------------------------------------------------------------------------

        Logger::logger->info("Making matrix A");
        auto A = makeA(_descriptors, dataToBeFit);
        Logger::logger->info("Done making matrix A");

        //Logger::logger->debug(format("{}|{}|{}|{}", A(0,0), A(0, 1), A(1,0), A(1, 1)));
        //Logger::logger->debug(format("A{}|{}|{}", A(0,0), A(1, 0), A(10,0)));
        //Logger::logger->debug(format("A{}|{}", A(0,0), A(1, 0)));
        // Logger::logger->info(matrixToString(A));
        // Logger::logger->debug(format("A{}", A.norm()));

        Logger::logger->info("Making feature vector b");
        auto b = makeB(_descriptors, dataToBeFit);
        Logger::logger->info("Done making feature vector b");

        //// ----------------------------------------------------------------------------------------------------

        Logger::logger->info("Doing linear algebra");

        Logger::logger->debug("Init Eigen::HouseholderQR");
        const Eigen::HouseholderQR<Eigen::MatrixXd> qr(A);

        Logger::logger->debug("Q^t");
        auto Qt = qr.householderQ().transpose();
        //Eigen::MatrixXd _Qtt = Qt;

        Logger::logger->debug("Q^t * b");
        auto Qt_b =  Qt * b;
         //Logger::logger->info(matrixToString(qr.householderQ().transpose()));
         //Logger::logger->info(matrixToString(_Qtt*_Qtt.transpose()));

        Logger::logger->debug("R");
        auto R = qr.matrixQR().topLeftCorner(A.cols(), A.cols());
        // Logger::logger->info(matrixToString(R));

        Logger::logger->debug("R^-1 * Qt_b");
        Eigen::VectorXd c = R.triangularView<Eigen::Upper>().solve(Qt_b.head(A.cols()));
        //Logger::logger->debug(format("{} {}", c[0], c[1]));

        Logger::logger->info("Finished linear algebra");

        size_t counter = 0;
        for (auto &descriptor: _descriptors) {
            size_t n = descriptor->nSparsePoints();

            Logger::logger->debug(format("desc AAAA{}", n));
            auto bb = vector<double>{c.data() + counter, c.data() + counter + n};
            Logger::logger->debug(format("desc BBBB{}", bb.size()));
            descriptor->setCoefficients(vector<double>{c.data() + counter, c.data() + counter + n});
            Logger::logger->debug(format("desc CCC{}", bb.size()));

            counter += n;
        }

        return make_shared<JgapPotential>(_descriptors);
    }

    vector<AtomicStructure> InRamJgapFit::subtractExternalContributions(const vector<AtomicStructure> &originalData) const {

        Logger::logger->info("Subtracting external potential contributions");

        vector dataToBeFit(originalData);

        // Subtract external potential contributions
        for (auto& structure: dataToBeFit) {
            for (const auto& potential: _externalPotentials) {
                const auto potentialsContribution = potential->predict(structure);
                structure.adjust(potentialsContribution, true, false);
            }
        }

        return dataToBeFit;
    }

    Eigen::MatrixXd InRamJgapFit::makeA(const vector<shared_ptr<Descriptor>> &descriptors,
                                        vector<AtomicStructure> &atomicStructures) const {

        size_t r = 0;
        vector<pair<size_t, AtomicStructure>> startingRowsK_nm;
        for (const auto& structure : atomicStructures) {
            startingRowsK_nm.emplace_back(r, structure);
            r += 1 + 3 * structure.atoms.size();
        }

        vector<tuple<size_t/*row*/, size_t/*col*/, size_t/*desc_idx*/>> startingPointsK_mm;
        size_t c = 0;
        for (size_t i = 0; i < descriptors.size(); i++) {
            startingPointsK_mm.emplace_back(r + c, c, i);
            c += descriptors[i]->nSparsePoints();
        }

        Logger::logger->info("Forming in-memory " + to_string(r+c) + "x" + to_string(c) + " A matrix");
        Eigen::MatrixXd resultingA(r + c, c);

        atomic counter(0);
        tbb::parallel_for_each(
            startingRowsK_nm.begin(), startingRowsK_nm.end(),
            [&](const pair<size_t, AtomicStructure>& structId) {

                size_t progress = ++counter;
                if (progress % (startingRowsK_nm.size() / 100) == 0) {
                    Logger::logger->debug(format(
                            "K_nm matrix formation progress: {} of {} ({}%)",
                            progress,
                            startingRowsK_nm.size(),
                            progress * 100 / startingRowsK_nm.size()
                            ));
                }

                fillInverseRootSigmaK_nm(descriptors, structId.second, resultingA, structId.first);
            }
        );

        tbb::parallel_for_each(
            startingPointsK_mm.begin(), startingPointsK_mm.end(),
            [&](const tuple<size_t/*row*/, size_t/*col*/, size_t/*desc_idx*/>& descriptorId) {
                Logger::logger->debug(format("K_mm for descriptor {}", get<2>(descriptorId)));
                fillU_mm(get<0>(descriptorId), get<1>(descriptorId), *descriptors[get<2>(descriptorId)], resultingA);
            }
        );

        return resultingA;
    }

    Eigen::VectorXd InRamJgapFit::makeB(const vector<shared_ptr<Descriptor>> &descriptors,
                                        const vector<AtomicStructure> &atomicStructures) {
        vector<double> b;
        for (auto& structure: atomicStructures) {
            if (structure.energy.has_value()) {
                b.push_back(structure.energy.value() * pow(structure.energySigma.value(), -1));
            }
            for (auto& atom: structure.atoms) {
                if (atom.force.has_value()) {
                    b.push_back(atom.force.value().x / atom.forceSigmas.value().x);
                    b.push_back(atom.force.value().y / atom.forceSigmas.value().y);
                    b.push_back(atom.force.value().z / atom.forceSigmas.value().z);
                }
            }
        }

        for (auto& descriptor: descriptors) {
            b.resize(b.size() + descriptor->nSparsePoints());
        }

        return Eigen::Map<Eigen::VectorXd>(b.data(), b.size());
    }

    void InRamJgapFit::fillInverseRootSigmaK_nm(const vector<shared_ptr<Descriptor>> &descriptors,
        const AtomicStructure &atomicStructure, Eigen::MatrixXd &A, size_t startingRow) {

        size_t currentColumn = 0;
        for (const auto& descriptor : descriptors) {
            auto contribution = descriptor->covariate(atomicStructure);

            for (size_t colInc = 0; colInc < contribution.size(); colInc++) {
                double contr = contribution[colInc].total; // debug
                A(startingRow, currentColumn + colInc) = contr * pow(atomicStructure.energySigma.value(), -1);
                //Logger::logger->debug(format("E_{}: {}", colInc,A(startingRow, currentColumn + colInc)));

                for (size_t rowInc = 0; rowInc < contribution[colInc].derivatives.size(); rowInc++) {
                    const auto derivative = contribution[colInc].derivatives[rowInc]; // FORCE IS NEGATIVE
                    const auto fSigmas = atomicStructure.atoms[rowInc].forceSigmas.value();
                    A(startingRow + rowInc * 3 + 1, currentColumn + colInc) = - derivative.x / fSigmas.x;
                    A(startingRow + rowInc * 3 + 2, currentColumn + colInc) = - derivative.y / fSigmas.y;
                    A(startingRow + rowInc * 3 + 3, currentColumn + colInc) = - derivative.z / fSigmas.z;

                    /*Logger::logger->debug(format("f_{}: {}|{}|{}", rowInc,
                    A(startingRow + rowInc * 3 + 1, currentColumn + colInc),
                    A(startingRow + rowInc * 3 + 2, currentColumn + colInc),
                    A(startingRow + rowInc * 3 + 3, currentColumn + colInc)
                    ));*/
                }
            }
            currentColumn += contribution.size();
        }
    }

    void InRamJgapFit::fillU_mm(size_t startingRow, size_t startingCol, Descriptor &descriptor,
        Eigen::MatrixXd &A) const {
        for (auto &[sparseId, K_mmPart]: descriptor.selfCovariate()) {

            size_t n = K_mmPart->rows();
            for (size_t i = 0; i < n; i++) (*K_mmPart)(i, i) += _jitter;

            auto K_mmConverted = convertToEigen(*K_mmPart);
            auto U_mmPart = choleskyDecomposition(K_mmConverted);

            for (size_t i = 0; i < n; i++) {
                for (size_t j = 0; j < n; j++) {
                    A(startingRow + sparseId + i, startingCol + sparseId + j) = U_mmPart(i, j);
                }
            }
        }
    }

    Eigen::MatrixXd InRamJgapFit::choleskyDecomposition(Eigen::MatrixXd &matrix) {
        Eigen::LLT<Eigen::MatrixXd> llt(matrix);
        if (llt.info() != Eigen::Success)
            Logger::logger->error("Cholesky decomposition failed: matrix not positive definite", true);

        return llt.matrixU();
    }

    Eigen::MatrixXd InRamJgapFit::convertToEigen(MatrixBlock &matrixBlock) {
        return Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>(
            matrixBlock.rawData().data(), matrixBlock.rows(), matrixBlock.columns()
            );
    }
}
