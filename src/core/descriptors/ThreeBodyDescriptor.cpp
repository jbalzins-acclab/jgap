#include "core/descriptors/ThreeBodyDescriptor.hpp"

#include <random>
#include <tbb/parallel_for_each.h>

#include "core/descriptors/kernels/ThreeBodySE.hpp"
#include "io/log/StdoutLogger.hpp"
#include "utils/Utils.hpp"

namespace jgap {
    ThreeBodyDescriptor::ThreeBodyDescriptor(const nlohmann::json &params) {
        if (params["kernel"]["cutoff"].is_number()) {
            _cutoff = params["kernel"]["cutoff"];
        } else {
            _cutoff = params["kernel"]["cutoff"]["cutoff"];
        }

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
                    _sparsePointsPerSpeciesTriplet[speciesTriplet].emplace_back(
                        triplet["x"], triplet["y"], triplet["z"]
                    );
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
            _sparsifier = ParserRegistry<PerSpecies3bSparsifier>::get(params["sparsify"]);
        }
    }

    nlohmann::json ThreeBodyDescriptor::serialize() {
        nlohmann::json sparseData;

        size_t counter = 0;
        for (auto &[speciesTriplet, sparseVectors]: _sparsePointsPerSpeciesTriplet) {
            vector<map<string, double>> converted;
            for (auto &vec: sparseVectors) {
                converted.push_back({
                        {"x", vec.x},
                        {"y", vec.y},
                        {"z", vec.z},
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

        return {
            {"sparse_data", sparseData},
            {"kernel", kernelData}
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

        _sparsePointsPerSpeciesTriplet = _sparsifier->sparsifyFromData(fromData);
        _coefficients.clear();
    }

    vector<Covariance> ThreeBodyDescriptor::covariate(const AtomicStructure &atomicStructure) {
        auto covariates = vector<Covariance>();

        auto indexMap = doIndex(atomicStructure);

        for (const auto& [speciesTriplet, sparsePoints]: _sparsePointsPerSpeciesTriplet) {
            for (const Vector3 sparsePoint : sparsePoints) {

                const auto u = _kernel->covariance(atomicStructure, indexMap[speciesTriplet], sparsePoint);
                const auto ddrs = _kernel->derivatives(atomicStructure, indexMap[speciesTriplet], sparsePoint);

                covariates.push_back({u, ddrs});
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

            auto tripletEnergies = vector(params.grid3b.size(), 0.0);

            vector<size_t> grid3bIndexes{};
            for (size_t i = 0; i < params.grid3b.size(); i++) {
                grid3bIndexes.push_back(i);
            }
            tbb::parallel_for_each(grid3bIndexes.begin(), grid3bIndexes.end(), [&](size_t iGrid) {

                const Vector3 invariantTriplet{
                    .x = params.grid3b[iGrid].x + params.grid3b[iGrid].y,
                    .y = pow(params.grid3b[iGrid].x - params.grid3b[iGrid].y, 2),
                    .z = sqrt(
                        pow(params.grid3b[iGrid].x, 2)
                        + pow(params.grid3b[iGrid].y, 2)
                        - 2.0 * params.grid3b[iGrid].x * params.grid3b[iGrid].x * params.grid3b[iGrid].z
                        )
                };

                for (size_t indexSparse = 0; indexSparse < sparsePoints.size(); indexSparse++) {
                    tripletEnergies[iGrid] += coefficients[indexSparse]
                                                * _kernel->covariance(invariantTriplet, sparsePoints[indexSparse]);
                }
            });

            result.tripletEnergies[speciesTriplet] = tripletEnergies;
        }

        return result;
    }

    map<SpeciesTriplet, ThreeBodyKernelIndex> ThreeBodyDescriptor::doIndex(
        const AtomicStructure &atomicStructure) const {

        map<SpeciesTriplet, ThreeBodyKernelIndex> indexes;

        for (size_t atomIndex = 0; atomIndex < atomicStructure.atoms.size(); atomIndex++) {
            auto atom = atomicStructure.atoms[atomIndex];

            // "nl" = neighbourList
            for (size_t nlIndex1 = 0; nlIndex1 < atom.neighbours->size(); nlIndex1++) {
                auto neighbour1 = atom.neighbours->at(nlIndex1);
                if (neighbour1.distance > _cutoff) continue;

                for (size_t nlIndex2 = nlIndex1 + 1; nlIndex2 < atom.neighbours->size(); nlIndex2++) {
                    auto neighbour2 = atom.neighbours->at(nlIndex2);
                    if (neighbour2.distance > _cutoff) continue;

                    auto speciesTriplet = SpeciesTriplet{atom.species,{
                        atomicStructure.atoms[neighbour1.index].species,
                        atomicStructure.atoms[neighbour2.index].species
                    }};

                    if (!indexes.contains(speciesTriplet)) {
                        indexes[speciesTriplet] = ThreeBodyKernelIndex();
                    }

                    indexes[speciesTriplet].push_back({atomIndex, nlIndex1, nlIndex2});
                }
            }
        }

        return indexes;

    }
}
