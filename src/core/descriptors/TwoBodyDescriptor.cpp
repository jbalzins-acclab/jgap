#include "core/descriptors/TwoBodyDescriptor.hpp"

#include <random>
#include <set>
#include <nlohmann/json.hpp>
#include <tbb/parallel_for_each.h>

#include "core/descriptors/kernels/TwoBodySE.hpp"
#include "io/log/StdoutLogger.hpp"
#include "io/parse/ParserRegistry.hpp"
#include "utils/Utils.hpp"

using namespace std;

namespace jgap {
    TwoBodyDescriptor::TwoBodyDescriptor(const nlohmann::json& params) {
        // Either explicitly in kernel or in cutoff obj specs
        if (params["kernel"]["cutoff"].is_number()) {
            _cutoff = params["kernel"]["cutoff"];
        } else {
            _cutoff = params["kernel"]["cutoff"]["cutoff"];
        }

        _kernel = ParserRegistry<TwoBodyKernel>::get(params["kernel"]);

        _sparsePointsPerSpeciesPair = {};
        _coefficients = {};
        if (params.contains("sparse_data")) {
            _sparsifier = nullptr;

            size_t nPts = 0;
            for (auto &[speciesStr, perPairData]: params["sparse_data"].items()) {
                SpeciesPair speciesPair = {split(speciesStr, ',')[0], split(speciesStr, ',')[1]};
                _sparsePointsPerSpeciesPair[speciesPair] = {};

                for (const auto &distance: perPairData["sparse_points"]) {
                    _sparsePointsPerSpeciesPair[speciesPair].push_back(distance);
                    nPts++;
                }

                if (perPairData.contains("coefficients")) {
                    for (const auto &coeff: perPairData["coefficients"]) {
                        _coefficients.push_back(coeff);
                    }
                }
            }

            if (!_coefficients.empty() && nPts != _coefficients.size()) {
                CurrentLogger::get()->error("Number of coefficients doesn't match number of 2b sparse points", true);
            }
        } else {
            _sparsifier = ParserRegistry<PerSpecies2bSparsifier>::get(params["sparsify"]);
        }
    }

    nlohmann::json TwoBodyDescriptor::serialize() {

        nlohmann::json sparseData;
        size_t counter = 0;
        for (auto &[speciesPair, sparsePoints] : _sparsePointsPerSpeciesPair) {
            sparseData[speciesPair.toString()] = {};
            sparseData[speciesPair.toString()]["sparse_points"] = sparsePoints;
            if (_coefficients.size() == nSparsePoints()) {
                sparseData[speciesPair.toString()]["coefficients"] = vector(
                    _coefficients.begin() + counter, _coefficients.begin() + counter + sparsePoints.size()
                );
            }
            counter += sparsePoints.size();
        }

        auto kernelData = _kernel->serialize();
        kernelData["type"] = _kernel->getType();

        return {
            {"sparse_data", sparseData},
            {"kernel", kernelData}
        };
    }

    size_t TwoBodyDescriptor::nSparsePoints() {
        size_t result = 0;
        for (auto &points: _sparsePointsPerSpeciesPair | views::values) {
            result += points.size();
        }
        return result;
    }

    void TwoBodyDescriptor::setSparsePoints(const vector<AtomicStructure> &fromData) {
        if (_sparsifier == nullptr) {
            CurrentLogger::get()->error("2b sparsifier not set", true);
        }

        _sparsePointsPerSpeciesPair = _sparsifier->sparsifyFromData(fromData);
        _coefficients.clear();
    }

    vector<Covariance> TwoBodyDescriptor::covariate(const AtomicStructure &atomicStructure) {
        auto covariates = vector<Covariance>();

        auto indexes = doIndex(atomicStructure);

        for (const auto& [speciesPair, sparsePoints]: _sparsePointsPerSpeciesPair) {
            for (const double sparsePoint : sparsePoints) {

                auto u = _kernel->covariance(atomicStructure, indexes[speciesPair], sparsePoint);
                auto ddrs = _kernel->derivatives(atomicStructure, indexes[speciesPair], sparsePoint);

                covariates.push_back({u, ddrs});
            }
        }

        return covariates;
    }

    vector<pair<size_t, shared_ptr<MatrixBlock>>> TwoBodyDescriptor::selfCovariate() {
        vector<pair<size_t, shared_ptr<MatrixBlock>>> result;

        size_t counter = 0;
        for (auto &sparsePoints: _sparsePointsPerSpeciesPair | views::values) {

            auto covariance = make_shared<MatrixBlock>(sparsePoints.size(), sparsePoints.size());

            for (size_t i = 0; i < sparsePoints.size(); i++) {
                for (size_t j = 0; j < sparsePoints.size(); j++) {
                    (*covariance)(i, j) = _kernel->covariance(sparsePoints[i], sparsePoints[j]); // TODO: cutoff??
                }
            }

            result.emplace_back(counter, covariance);
            counter += sparsePoints.size();
        }

        return result;
    }

    TabulationData TwoBodyDescriptor::tabulate(const TabulationParams &params) {

        TabulationData result{};

        size_t counter = 0;
        for (const auto &[speciesPair, sparsePoints]: _sparsePointsPerSpeciesPair) {
            vector coefficients(_coefficients.begin()+counter, _coefficients.begin()+counter + sparsePoints.size());
            counter += sparsePoints.size();

            auto pairEnergies = vector(params.grid2b.size(), 0.0);

            for (size_t iGrid = 0; iGrid < params.grid2b.size(); iGrid++) {
                for (size_t indexSparse = 0; indexSparse < sparsePoints.size(); indexSparse++) {
                    pairEnergies[iGrid] += coefficients[indexSparse]
                                        * _kernel->covariance(params.grid2b[iGrid], sparsePoints[indexSparse]);
                }
            }

            result.pairEnergies[speciesPair] = pairEnergies;
        }

        return result;
    }

    map<SpeciesPair, TwoBodyKernelIndex> TwoBodyDescriptor::doIndex(const AtomicStructure &atomicStructure) const {

        map<SpeciesPair, TwoBodyKernelIndex> indexes;

        for (size_t atomIndex = 0; atomIndex < atomicStructure.atoms.size(); atomIndex++) {
            auto atom = atomicStructure.atoms[atomIndex];

            for (size_t neighbourListIndex = 0; neighbourListIndex < atom.neighbours->size(); neighbourListIndex++) {
                auto neighbour = atom.neighbours->at(neighbourListIndex);

                if (neighbour.index < atomIndex) continue;
                if (neighbour.distance > _cutoff) continue;

                auto speciesPair = SpeciesPair{atom.species, atomicStructure.atoms[neighbour.index].species};
                if (!indexes.contains(speciesPair)) {
                    indexes[speciesPair] = TwoBodyKernelIndex();
                }

                indexes[speciesPair].push_back({atomIndex, neighbourListIndex});
            }
        }

        return indexes;
    }
}
