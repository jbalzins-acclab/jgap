#include "core/fit/InRamJgapFit.hpp"

#include "core/matrices/sigmas/SimpleSigmaRules.hpp"
#include "core/neighbours/NeighbourFinder.hpp"
#include "io/log/StdoutLogger.hpp"
#include "utils/Utils.hpp"

#include <tbb/parallel_for_each.h>

namespace jgap {

    InRamJgapFit::InRamJgapFit(const nlohmann::json &params) {

        _descriptors = {};
        for (const auto& [label, descriptor] : params["descriptors"].items()) {
            _descriptors[label] = (ParserRegistry<Descriptor>::get(descriptor));
        }

        _sigmaRules = ParserRegistry<SigmaRules>::get(params["sigma_rules"]);
        _jitter = params.value("jitter", 1e-8);
    }

    shared_ptr<Potential> InRamJgapFit::fit(const vector<AtomicStructure>& trainingData) {
        CurrentLogger::get()->info("Starting JGAP fit");

        vector _trainingData(trainingData);

        double maxCutoff = 0;

        CurrentLogger::get()->info("Checking sparse points");

        auto descriptorsAsVec = vector<shared_ptr<Descriptor>>();
        for (const auto& descriptor: _descriptors | views::values) {
            descriptorsAsVec.push_back(descriptor);

            // To avoid ugly "cutoff" in sparse json
            NeighbourFinder::findNeighbours(_trainingData, descriptor->getCutoff());
            if (descriptor->nSparsePoints() == 0) {
                CurrentLogger::get()->info(format(
                    "No sparse points found in descriptor {} -> setting from data",
                    descriptor->serialize().dump()
                    ));
                descriptor->setSparsePoints(_trainingData);
            }

            maxCutoff = max(maxCutoff, descriptor->getCutoff());
        }
        CurrentLogger::get()->info("Sparsification complete");

        CurrentLogger::get()->debug("Full neighbour-list");
        NeighbourFinder::findNeighbours(_trainingData, maxCutoff);

        CurrentLogger::get()->info("Per-structure regularization setup");
        for (auto& structure: _trainingData) {
            _sigmaRules->fillSigmas(structure);
        }
        CurrentLogger::get()->info("Per-structure regularization setup finished");

        //// ----------------------------------------------------------------------------------------------------

        CurrentLogger::get()->info("Making matrix A");
        auto A = makeA(descriptorsAsVec, _trainingData);
        CurrentLogger::get()->info("Done making matrix A");
        // CurrentLogger::get()->info(matrixToString(A));

        CurrentLogger::get()->info("Making feature vector b");
        auto b = makeB(descriptorsAsVec, _trainingData);
        CurrentLogger::get()->info("Done making feature vector b");

        //// ----------------------------------------------------------------------------------------------------

        CurrentLogger::get()->info("Doing linear algebra");

        CurrentLogger::get()->debug("Init Eigen::HouseholderQR");
        const Eigen::HouseholderQR<Eigen::MatrixXd> qr(A);

        CurrentLogger::get()->debug("Q^t");
        auto Qt = qr.householderQ().transpose();

        CurrentLogger::get()->debug("Q^t * b");
        auto Qt_b =  Qt * b;

        CurrentLogger::get()->debug("R");
        auto R = qr.matrixQR().topLeftCorner(A.cols(), A.cols());

        CurrentLogger::get()->debug("R^-1 * Qt_b");
        Eigen::VectorXd c = R.triangularView<Eigen::Upper>().solve(Qt_b.head(A.cols()));
        CurrentLogger::get()->debug("c" + to_string(c[0]));

        CurrentLogger::get()->info("Finished linear algebra");

        size_t counter = 0;
        for (const auto &descriptor: _descriptors | views::values) {
            size_t n = descriptor->nSparsePoints();

            // auto bb = vector<double>{c.data() + counter, c.data() + counter + n};
            descriptor->setCoefficients(vector<double>{c.data() + counter, c.data() + counter + n});

            counter += n;
        }

        return make_shared<JgapPotential>(_descriptors);
    }

    Eigen::MatrixXd InRamJgapFit::makeA(const vector<shared_ptr<Descriptor>> &descriptors,
                                        const vector<AtomicStructure> &atomicStructures) const {

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

        CurrentLogger::get()->info("Forming in-memory " + to_string(r+c) + "x" + to_string(c) + " A matrix");
        Eigen::MatrixXd resultingA(r + c, c);

        atomic counter(0);
        tbb::parallel_for_each(
            startingRowsK_nm.begin(), startingRowsK_nm.end(),
            [&](const pair<size_t, AtomicStructure>& structId) {

                size_t progress = ++counter;
                if (progress % max(startingRowsK_nm.size() / 100, 1uz) == 0) {
                    CurrentLogger::get()->debug(format(
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
                CurrentLogger::get()->debug(format("K_mm for descriptor {}", get<2>(descriptorId)));
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
                //CurrentLogger::get()->debug(format("E_{}: {}", colInc,A(startingRow, currentColumn + colInc)));
                for (size_t rowInc = 0; rowInc < contribution[colInc].derivatives.size(); rowInc++) {
                    const auto derivative = contribution[colInc].derivatives[rowInc]; // FORCE IS NEGATIVE
                    const auto fSigmas = atomicStructure.atoms[rowInc].forceSigmas.value();
                    A(startingRow + rowInc * 3 + 1, currentColumn + colInc) = - derivative.x / fSigmas.x;
                    A(startingRow + rowInc * 3 + 2, currentColumn + colInc) = - derivative.y / fSigmas.y;
                    A(startingRow + rowInc * 3 + 3, currentColumn + colInc) = - derivative.z / fSigmas.z;
                CurrentLogger::get()->debug(format("U_ffX{}:{}:{}", derivative.x, startingRow + rowInc * 3 + 1,  currentColumn + colInc));
                CurrentLogger::get()->debug(format("U_ffY{}:{}:{}", derivative.y, startingRow + rowInc * 3 + 2,  currentColumn + colInc));
                CurrentLogger::get()->debug(format("U_ffZ{}:{}:{}", derivative.z, startingRow + rowInc * 3 + 3,  currentColumn + colInc));

                    /*CurrentLogger::get()->debug(format("f_{}: {}|{}|{}", rowInc,
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
            CurrentLogger::get()->error("Cholesky decomposition failed: matrix not positive definite", true);

        return llt.matrixU();
    }

    Eigen::MatrixXd InRamJgapFit::convertToEigen(MatrixBlock &matrixBlock) {
        return Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>(
            matrixBlock.rawData().data(), matrixBlock.rows(), matrixBlock.columns()
            );
    }
}
