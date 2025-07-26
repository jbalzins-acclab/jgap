//
// Created by Jegors Balzins on 18.6.2025.
//

#include "core/descriptors/EamDescriptor.hpp"

#include <random>

#include "core/descriptors/kernels/EamSE.hpp"
#include "io/log/StdoutLogger.hpp"
#include "io/parse/ParserRegistry.hpp"

namespace jgap {
    EamDescriptor::EamDescriptor(const nlohmann::json &params) {
        CurrentLogger::get()->debug("Parsing EAM descriptor params");

        _sparsePointsPerSpecies = {};
        _coefficients = {};

        _kernel = ParserRegistry<EamKernel>::get(params["kernel"]);

        if (params.contains("sparse_data")) {
            _sparsifier = nullptr; // TODO ?

            size_t nPts = 0;
            for (const auto &[species, sparseData] : params["sparse_data"].items()) {
                _sparsePointsPerSpecies[species] = {};
                for (const auto& density: sparseData["sparse_points"]) {
                    _sparsePointsPerSpecies[species].push_back(density);
                    nPts++;
                }

                if (sparseData.contains("coefficients")) {
                    for (const auto& coeff: sparseData["coefficients"]) {
                        _coefficients.push_back(coeff);
                    }
                }
            }

            if (!_coefficients.empty() && nPts != _coefficients.size()) {
                CurrentLogger::get()->error("Coefficients can't be provided partially for EAM descriptor", true);
            }
        } else {
            _sparsifier = ParserRegistry<Sparsifier>::get(params["sparsify"]);
        }

        _maxCutoff = 0;
        for (const auto& pfParams: params["pair_functions"]) {

            auto pf = ParserRegistry<EamPairFunction>::get(pfParams);

            // TODO: error on specified double
            if (pfParams.contains("species")) {
                auto s1 = pfParams["species"][0];
                auto s2 = pfParams["species"][1];
                _pairFunctions[{s1, s2}] = pf;
                _pairFunctions[{s2, s1}] = pf;
            } else if (pfParams.contains("species_ordered")) {
                auto s1 = pfParams["species_ordered"][0];
                auto s2 = pfParams["species_ordered"][1];
                _pairFunctions[{s1, s2}] = pf;
            } else { // default
                _defaultPairFunction = pf;
            }

            _maxCutoff = max(_maxCutoff, pfParams["cutoff"].get<double>());
        }
    }

    nlohmann::json EamDescriptor::serialize() {

        nlohmann::json sparseData{};
        size_t counter = 0;

        for (const auto &[species, sparseDensities] : _sparsePointsPerSpecies) {
            sparseData[species] = {
                {"sparse_points", sparseDensities},
                {"coefficients", vector(
                    _coefficients.begin() + counter,
                    _coefficients.begin() + counter + sparseDensities.size()
                    )
                }
            };
            counter += sparseDensities.size();
        }

        auto kernelData = _kernel->serialize();
        kernelData["type"] = _kernel->getType();

        nlohmann::json pfData = nlohmann::json::array();

        auto defaultPairFunctionData = _defaultPairFunction->serialize();
        defaultPairFunctionData["type"] = _defaultPairFunction->getType();
        pfData.push_back(defaultPairFunctionData);

        for (const auto& [orderedSpeciesPair, pf]: _pairFunctions) {
            auto newPfData = pf->serialize();
            newPfData["type"] = pf->getType();
            newPfData["species_ordered"] = vector{orderedSpeciesPair.first, orderedSpeciesPair.second};
            pfData.push_back(newPfData);
        }

        return {
            {"sparse_data", sparseData},
            {"kernel", kernelData},
            {"pair_functions", pfData}
        };
    }

    void EamDescriptor::setSparsePoints(const vector<AtomicStructure> &fromData) {
        if (_sparsifier == nullptr) {
            CurrentLogger::get()->error("EAM sparsifier not set", true);
        }
        CurrentLogger::get()->info("Doing EAM sparsification from data");

        map<Species, vector<vector<double>>> allDensitiesPerSpecies;
        vector<EamKernelIndex> indexArr;
        for (const auto& structure : fromData) {
            for (auto structureIndex = doIndex(structure);
                 const auto& [species, densities]: structureIndex) {
                if (!allDensitiesPerSpecies.contains(species)) {
                    allDensitiesPerSpecies[species] = {};
                }
                for (const auto& densityData: densities) {
                    if (densityData.density > 3) {
                        allDensitiesPerSpecies[species].push_back(vector{densityData.density});
                    }
                }
            }
        }

        _sparsePointsPerSpecies = {};
        for (const auto& [species, densities]: allDensitiesPerSpecies) {
            _sparsePointsPerSpecies[species] = {};
            for (const vector<double>& density: _sparsifier->selectSparsePoints(densities)) {
                _sparsePointsPerSpecies[species].push_back(density[0]);
            }
        }
        _coefficients.clear();
    }

    size_t EamDescriptor::nSparsePoints() {
        size_t result = 0;
        for (const auto &densities: _sparsePointsPerSpecies | views::values) {
            result += densities.size();
        }
        return result;
    }

    vector<Covariance> EamDescriptor::covariate(const AtomicStructure &atomicStructure) {
        vector<Covariance> result;

        EamKernelIndex kernelIndex = doIndex(atomicStructure);

        for (auto &[species, sparseDensities]: _sparsePointsPerSpecies) {
            if (!kernelIndex.contains(species)) kernelIndex[species] = {};

            for (double sparseDensity: sparseDensities) {
                result.push_back(
                    _kernel->covariance(atomicStructure, kernelIndex.at(species), sparseDensity)
                    );
            }
        }

        return result;
    }

    vector<pair<size_t, shared_ptr<MatrixBlock>>> EamDescriptor::selfCovariate() {

        vector<pair<size_t, shared_ptr<MatrixBlock>>> result;
        size_t startingRC = 0;

        for (const auto &densities: _sparsePointsPerSpecies | views::values) {

            auto elementBlock = make_shared<MatrixBlock>(densities.size(), densities.size());
            for (size_t i = 0; i < densities.size(); i++) {
                for (size_t j = 0; j < densities.size(); j++) {
                    (*elementBlock)(i, j) = _kernel->covariance(densities[i], densities[j]);
                }
            }

            result.emplace_back(startingRC, elementBlock);
            startingRC += densities.size();
        }

        return result;
    }

    TabulationData EamDescriptor::tabulate(const TabulationParams &params) {
        EamTabulationData result;

        result.rhoMax = 0.0;
        for (const auto& points: _sparsePointsPerSpecies | views::values) {
            result.rhoMax = max(result.rhoMax, ranges::max(points));
        }
        result.rhoMax += 3.5; // TODO: this is very SqExp specific

        const double rhoStep = result.rhoMax / static_cast<double>(params.nDensities-1);

        size_t counter = 0;
        for (const auto& [species, points]: _sparsePointsPerSpecies) {
            vector energiesPerSpecies(params.nDensities, 0.0);
            vector coefficients(_coefficients.begin()+counter, _coefficients.begin()+counter + points.size());
            counter += points.size();

            for (size_t iGrid = 0; iGrid < params.nDensities; iGrid++) {
                double density = rhoStep * static_cast<double>(iGrid);

                for (size_t indexSparse = 0; indexSparse < points.size(); indexSparse++) {
                    energiesPerSpecies[iGrid] += coefficients[indexSparse]
                                                  * _kernel->covariance(points[indexSparse], density);
                }
            }

            result.embeddingEnergies[species] = energiesPerSpecies;
        }

        for (const auto& species1: _sparsePointsPerSpecies | views::keys) {
            for (const auto& species2: _sparsePointsPerSpecies | views::keys) {

                auto speciesPair = OrderedSpeciesPair{species1, species2};

                auto pairFunction = _defaultPairFunction;
                if (_pairFunctions.contains(speciesPair)) {
                    pairFunction = _pairFunctions[speciesPair];
                }

                vector<double> pairDensities{};
                for (const double& gridDensity: params.grid2b) {
                    pairDensities.push_back(pairFunction->evaluate(gridDensity));
                }

                result.eamDensities[speciesPair] = pairDensities;
            }
        }

        TabulationData resultFull{};
        resultFull.eamTabulationData = vector{result}; // TODO? c++20
        return resultFull;
    }

    EamKernelIndex EamDescriptor::doIndex(const AtomicStructure &structure) const {

        EamKernelIndex result{};

        for (size_t atomIdx = 0; atomIdx < structure.size(); atomIdx++) {

            double totalDensity = 0;
            vector<pair<NeighbourData, double>> densityDerivatives;

            Species species = structure.species[atomIdx];
            for (NeighbourData neighbour: structure.neighbours.value()[atomIdx]) {
                if (neighbour.distance > _maxCutoff) continue;

                pair orderedSpeciesPair = {
                    structure.species[neighbour.index],  species
                };

                shared_ptr<EamPairFunction> pf;
                if (_pairFunctions.contains(orderedSpeciesPair)) {
                    pf = _pairFunctions.at(orderedSpeciesPair);
                } else if (_defaultPairFunction != nullptr) {
                    pf = _defaultPairFunction;
                } else {
                    continue;
                }

                totalDensity += pf->evaluate(neighbour.distance);
                densityDerivatives.emplace_back(
                    neighbour,
                    pf->differentiate(neighbour.distance)
                );
            }

            if (!result.contains( species)) {
                result[species] = {};
            }

            result[species].push_back({atomIdx, totalDensity, densityDerivatives});
        }

        return result;
    }
}
