#include "core/descriptors/ThreeBodyDescriptor.hpp"

#include <random>
#include <tbb/parallel_for_each.h>

#include "core/descriptors/kernels/ThreeBodySE.hpp"
#include "io/log/StdoutLogger.hpp"
#include "utils/Utils.hpp"

namespace jgap {
    ThreeBodyDescriptor::ThreeBodyDescriptor(const nlohmann::json &params) {
        CurrentLogger::get()->debug("Parsing 3b descriptor params");

        _cutoffFunction = ParserRegistry<CutoffFunction>::get(params["cutoff"]);
        _cutoff = _cutoffFunction->getCutoff();

        _kernel = ParserRegistry<ThreeBodyKernel>::get(params["kernel"]);

        _sparsePointsPerSpeciesTriplet = {};
        _coefficients = {};
        if (params.contains("sparse_data")) {
            _sparsifier = nullptr;

            size_t nPts = 0;
            for (auto &[speciesStr, perTripletData]: params["sparse_data"].items()) {
                SpeciesTriplet speciesTriplet = {
                    .root=split(speciesStr, ',')[0],
                    .nodes={split(speciesStr, ',')[1], split(speciesStr, ',')[2]}
                };
                _sparsePointsPerSpeciesTriplet[speciesTriplet] = {};

                for (const auto &triplet: perTripletData["sparse_points"]) {
                    Vector3 q = {triplet["x"], triplet["y"], triplet["z"]};
                    _sparsePointsPerSpeciesTriplet[speciesTriplet].emplace_back(q, invariantTripletToCutoff(q));
                    nPts++;
                }

                if (perTripletData.contains("coefficients")) {
                    for (const auto &coeff: perTripletData["coefficients"]) {
                        _coefficients.push_back(coeff);
                    }
                }
            }

            if (!_coefficients.empty() && nPts != _coefficients.size()) {
                CurrentLogger::get()->error("Number of coefficients doesn't match number of 3b sparse points", true);
            }

        } else {
            _sparsifier = ParserRegistry<Sparsifier>::get(params["sparsify"]);
        }
    }

    nlohmann::json ThreeBodyDescriptor::serialize() {
        nlohmann::json sparseData;

        size_t counter = 0;
        for (auto &[speciesTriplet, sparseVectors]: _sparsePointsPerSpeciesTriplet) {
            vector<map<string, double>> converted;
            for (auto &sparsePoint: sparseVectors) {
                converted.push_back({
                        {"x", sparsePoint.q.x},
                        {"y", sparsePoint.q.y},
                        {"z", sparsePoint.q.z},
                    });
            }
            sparseData[speciesTriplet.toString()] = {};
            sparseData[speciesTriplet.toString()]["sparse_points"] = converted;
            if (_coefficients.size() == nSparsePoints()) {
                sparseData[speciesTriplet.toString()]["coefficients"] = vector(
                    _coefficients.begin() + counter, _coefficients.begin() + counter + sparseVectors.size()
                );
            }
            counter += sparseVectors.size();
        }

        auto kernelData = _kernel->serialize();
        kernelData["type"] = _kernel->getType();

        auto cutoffData = _cutoffFunction->serialize();
        cutoffData["type"] = _cutoffFunction->getType();

        return {
            {"sparse_data", sparseData},
            {"kernel", kernelData},
            {"cutoff", cutoffData}
        };
    }

    size_t ThreeBodyDescriptor::nSparsePoints() {
        size_t result = 0;

        for (auto &sparseVectors: _sparsePointsPerSpeciesTriplet | views::values) {
            result += sparseVectors.size();
        }

        return result;
    }

    void ThreeBodyDescriptor::setSparsePoints(const vector<AtomicStructure> &fromData) {
        if (_sparsifier == nullptr) {
            CurrentLogger::get()->error("3b sparsifier not set", true);
        }
        CurrentLogger::get()->info("Doing 3b sparsification from data");

        map<SpeciesTriplet, vector<vector<double>>> allTriplets;
        for (const auto &structure: fromData) {
            for (const auto structureTriplets = doIndex(structure);
                 const auto &[speciesTriplet, points]: structureTriplets) {
                if (!allTriplets.contains(speciesTriplet)) {
                    allTriplets[speciesTriplet] = {};
                }

                for (const auto &point: points) {
                    allTriplets[speciesTriplet].push_back(vector{point.q.x, point.q.y, point.q.z});
                }
            }
        }

        _sparsePointsPerSpeciesTriplet.clear();
        _coefficients.clear();
        for (const auto &[speciesTriplet, allPoints]: allTriplets) {
            _sparsePointsPerSpeciesTriplet[speciesTriplet] = {};

            for (const vector<double> &point: _sparsifier->selectSparsePoints(allPoints)) {
                auto q = Vector3{point[0], point[1], point[2]};
                _sparsePointsPerSpeciesTriplet[speciesTriplet].push_back({
                    .q = q,
                    .fCut = invariantTripletToCutoff(q)
                });
            }
        }
    }

    void ThreeBodyDescriptor::setSparsePoints(
            const map<SpeciesTriplet, vector<Vector3>> &sparsePointsPerSpeciesTriplet) {
        _sparsePointsPerSpeciesTriplet.clear();

        for (const auto &[speciesTriplet, sparsePoints]: sparsePointsPerSpeciesTriplet) {
            _sparsePointsPerSpeciesTriplet[speciesTriplet] = {};

            for (const auto &point: sparsePoints) {
                auto q = Vector3{point.x, point.y, point.z};
                _sparsePointsPerSpeciesTriplet[speciesTriplet].push_back({
                    .q = q,
                    .fCut = invariantTripletToCutoff(q)
                });
            }
        }
    }

    vector<Covariance> ThreeBodyDescriptor::covariate(const AtomicStructure &atomicStructure) {
        auto covariates = vector<Covariance>();

        auto indexMap = doIndex(atomicStructure);

        for (const auto& [speciesTriplet, sparsePoints]: _sparsePointsPerSpeciesTriplet) {
            for (const auto& sparsePoint : sparsePoints) {
                covariates.push_back(_kernel->covariance(atomicStructure, indexMap[speciesTriplet], sparsePoint));
            }
        }

        return covariates;
    }

    vector<pair<size_t, shared_ptr<MatrixBlock>>> ThreeBodyDescriptor::selfCovariate() {
        vector<pair<size_t, shared_ptr<MatrixBlock>>> result;

        size_t counter = 0;
        for (auto &sparsePoints: _sparsePointsPerSpeciesTriplet | views::values) {

            auto covariance = make_shared<MatrixBlock>(sparsePoints.size(), sparsePoints.size());

            for (size_t i = 0; i < sparsePoints.size(); i++) {
                for (size_t j = 0; j < sparsePoints.size(); j++) {
                    (*covariance)(i, j) = _kernel->covariance(sparsePoints[i], sparsePoints[j]);
                }
            }

            result.emplace_back(counter, covariance);
            counter += sparsePoints.size();
        }

        return result;
    }

    TabulationData ThreeBodyDescriptor::tabulate(const TabulationParams &params) {

        TabulationData result{};

        size_t counter = 0;
        for (const auto &[speciesTriplet, sparsePoints]: _sparsePointsPerSpeciesTriplet) {
            vector coefficients(_coefficients.begin()+counter, _coefficients.begin()+counter + sparsePoints.size());
            counter += sparsePoints.size();

            auto tripletEnergies = vector(
                params.grid3b.size(),
                vector(params.grid3b[0].size(), vector(params.grid3b[0][0].size(), 0.0))
                );

            vector<array<size_t, 3>> grid3bIndexes{};
            for (size_t i = 0; i < params.grid3b.size(); i++) {
                for (size_t j = i; j < params.grid3b[i].size(); j++) {
                    for (size_t k = 0; k < params.grid3b[i][j].size(); k++) {
                        grid3bIndexes.push_back({i, j, k});
                    }
                }
            }
            tbb::parallel_for_each(grid3bIndexes.begin(), grid3bIndexes.end(), [&](const array<size_t, 3> &iGrid) {
                const Vector3 gridPoint = params.grid3b[iGrid[0]][iGrid[1]][iGrid[2]];

                const Vector3 invariantTriplet = toInvariantTriplet(
                    gridPoint.x,
                    gridPoint.y,
                    sqrt(max/*numeric safety*/(
                        pow(gridPoint.x, 2) + pow(gridPoint.y, 2) - 2.0 * gridPoint.x * gridPoint.y * gridPoint.z, 0.0
                        ))
                );

                for (size_t indexSparse = 0; indexSparse < sparsePoints.size(); indexSparse++) {
                    const double contribution = coefficients[indexSparse] * _kernel->covariance(
                        {.q = invariantTriplet, .fCut = invariantTripletToCutoff(invariantTriplet)},
                        sparsePoints[indexSparse]
                        ) * 2.0/*q_ijk + q_jik*/;
                    tripletEnergies[iGrid[0]][iGrid[1]][iGrid[2]] += contribution;
                }
                tripletEnergies[iGrid[1]][iGrid[0]][iGrid[2]] = tripletEnergies[iGrid[0]][iGrid[1]][iGrid[2]];
            });

            result.tripletEnergies[speciesTriplet] = tripletEnergies;
        }

        return result;
    }

    double ThreeBodyDescriptor::invariantTripletToCutoff(const Vector3 &t) const {
        const double dDiff = sqrt(t.y);
        const double d1 = (dDiff + t.x) / 2.0;
        const double d2 = (t.x - dDiff) / 2.0;

        return _cutoffFunction->evaluate(d1) * _cutoffFunction->evaluate(d2);
    }

    Vector3 ThreeBodyDescriptor::toInvariantTriplet(double r01, double r02, double r12) {
        return {r01 + r02, (r01-r02) * (r01 - r02), r12};
    }

    array<Vector3, 3> ThreeBodyDescriptor::invariantTripletGradients(double r01, double r02) {
        return {
            Vector3{1, 2 * (r01 - r02), 0},
            Vector3{1, 2 * (r02 - r01), 0},
            Vector3{0, 0, 1},
        };
    }

    map<SpeciesTriplet, ThreeBodyKernelIndex> ThreeBodyDescriptor::doIndex(
                                                const AtomicStructure &atomicStructure) const {

        map<SpeciesTriplet, ThreeBodyKernelIndex> indexes;

        for (size_t atomIndex = 0; atomIndex < atomicStructure.size(); atomIndex++) {
            auto atom0 = atomicStructure[atomIndex];

            // "nl" = neighbourList
            for (size_t nlIndex1 = 0; nlIndex1 < atom0.neighbours().size(); nlIndex1++) {
                auto neighbour1 = atom0.neighbours()[nlIndex1];
                auto atom1 = atomicStructure[neighbour1.index];
                if (neighbour1.distance > _cutoff) continue;

                for (size_t nlIndex2 = nlIndex1 + 1; nlIndex2 < atom0.neighbours().size(); nlIndex2++) {
                    auto neighbour2 = atom0.neighbours()[nlIndex2];
                    auto atom2 = atomicStructure[neighbour2.index];
                    if (neighbour2.distance > _cutoff) continue;

                    auto speciesTriplet = SpeciesTriplet{atom0.species(),{atom1.species(), atom2.species()}};

                    if (!indexes.contains(speciesTriplet)) {
                        indexes[speciesTriplet] = ThreeBodyKernelIndex();
                    }

                    array r_ij = {
                        atom1.position() + neighbour1.offset - atom0.position(),
                        atom2.position() + neighbour2.offset - atom0.position(),
                        atom2.position() + neighbour2.offset - (atom1.position() + neighbour1.offset)
                    };
                    indexes[speciesTriplet].push_back({
                        .atomIndex = {atomIndex, neighbour1.index, neighbour2.index},
                        .r_ij = r_ij,
                        .grad_rij_wrt_rj = {r_ij[0].normalize(), r_ij[1].normalize(), r_ij[2].normalize()},
                        .fCut01 = _cutoffFunction->evaluate(r_ij[0].len()),
                        .fCut02 = _cutoffFunction->evaluate(r_ij[1].len()),
                        .dfCut_dr_01 = _cutoffFunction->differentiate(r_ij[0].len()),
                        .dfCut_dr_02 = _cutoffFunction->differentiate(r_ij[1].len()),
                        .q = toInvariantTriplet(r_ij[0].len(), r_ij[1].len(), r_ij[2].len()),
                        .dq_k_dr_ij = invariantTripletGradients(r_ij[0].len(), r_ij[1].len())
                    });
                }
            }
        }

        return indexes;
    }
}
