#include "core/descriptors/eam/EamDescriptor.hpp"

#include <random>
#include <utility>

#include "core/descriptors/eam/EamSE.hpp"
#include "io/log/StdoutLogger.hpp"
#include "io/parse/ParserRegistry.hpp"

namespace jgap {
    EamDescriptor::EamDescriptor(vector<shared_ptr<EamKernel>> kernels,
                                 shared_ptr<EamPairFunction> defaultPairFunction,
                                 map<ContributorReceiverSpecies, shared_ptr<EamPairFunction>> pairFunctions)
        : _kernels(std::move(kernels)),
          _defaultPairFunction(std::move(defaultPairFunction)),
          _pairFunctions(std::move(pairFunctions)) {

        _maxCutoff = 0;
        for (const auto& pf: pairFunctions | views::values) {
            _maxCutoff = max(_maxCutoff, pf->getCutoff());
        }
        mapKernelIds();
    }

    EamDescriptor::EamDescriptor(const nlohmann::json &params) {
        CurrentLogger::get()->debug("Parsing EAM descriptor params");

        _kernels = {};
        if (params.contains("kernels")) {
            for (const nlohmann::json& kernelParams : params["kernels"]) {
                _kernels.push_back(ParserRegistry<EamKernel>::get(kernelParams));
            }
        }
        mapKernelIds();

        _kernelSetups = {};
        if (params.contains("kernel_setups")) {
            for (const nlohmann::json& kernelSetup : params["kernel_setups"]) {
                _kernelSetups.push(kernelSetup);
            }
        }
        if (params.contains("kernel_setup")) {
            _kernelSetups.push(params["kernel_setup"]);
        }

        _maxCutoff = 0;
        for (const auto& pfParams: params["pair_functions"]) {

            auto pf = ParserRegistry<EamPairFunction>::get(pfParams);

            if (pfParams.contains("species")) {
                auto s1 = pfParams["species"][0];
                auto s2 = pfParams["species"][1];
                _pairFunctions[{s1, s2}] = pf;
                _pairFunctions[{s2, s1}] = pf;
            } else if (pfParams.contains("contributor_receiver_species")) {
                auto s1 = pfParams["contributor_receiver_species"][0];
                auto s2 = pfParams["contributor_receiver_species"][1];
                _pairFunctions[{s1, s2}] = pf;
            } else {
                _defaultPairFunction = pf;
            }

            _maxCutoff = max(_maxCutoff, pf->getCutoff());
        }
    }

    nlohmann::json EamDescriptor::serialize() {

        auto kernelsData = nlohmann::json::array();
        for (const auto &kernel : _kernels) {
            kernelsData.push_back(kernel->serialize());
            kernelsData.back()["type"] = kernel->getType();
        }

        nlohmann::json pfData = nlohmann::json::array();

        auto defaultPairFunctionData = _defaultPairFunction->serialize();
        defaultPairFunctionData["type"] = _defaultPairFunction->getType();
        pfData.push_back(defaultPairFunctionData);

        for (const auto& [speciesPair, pf]: _pairFunctions) {
            auto newPfData = pf->serialize();
            newPfData["type"] = pf->getType();
            newPfData["contributor_receiver_species"] = vector{speciesPair.contributor, speciesPair.receiver};
            pfData.push_back(newPfData);
        }

        return {
            {"kernels", kernelsData},
            {"pair_functions", pfData}
        };
    }

    CutoffRanges EamDescriptor::getCutoff() {
        auto res = CutoffRanges{.twoBody = _maxCutoff};

        double minDensity = 0.0, maxDensity = numeric_limits<double>::min();
        for (const auto& kernel: _kernels) {
            minDensity = min(minDensity, kernel->getDensity());
            maxDensity = max(maxDensity, kernel->getDensity());
        }
        res.minEam = minDensity;
        res.maxEam = maxDensity;

        return res;
    }

    vector<shared_ptr<IKernel>> EamDescriptor::getKernels() {
        vector<shared_ptr<IKernel>> res;
        for (const auto& kernelIds : _kernelIndicesPerSpecies | views::values) {
            for (const auto& id: kernelIds) {
                res.push_back(_kernels[id]);
            }
        }
        return res;
    }

    void EamDescriptor::setupKernels(const vector<AtomicStructure> &fromData) {
        CurrentLogger::get()->info("Doing EAM sparsification from data");

        if (_kernelSetups.empty()) {
            CurrentLogger::get()->warn("All 3b kernels were pre-set");
            return;
        }

        map<Species, vector<vector<double>>> allDensitiesPerSpecies;
        vector<EamKernelIndex> indexArr;
        for (const auto& structure : fromData) {
            for (auto structureIndex = doIndex(structure);
                 const auto& [species, densities]: structureIndex) {
                if (!allDensitiesPerSpecies.contains(species)) {
                    allDensitiesPerSpecies[species] = {};
                }
                for (const auto& densityData: densities) {
                    // TODO: low density cap - avoid zero of isolated_atom density(?):
                    // if (densityData.density < 0.2) continue;
                    allDensitiesPerSpecies[species].push_back(vector{densityData.density});
                }
            }
        }

        while (!_kernelSetups.empty()) {
            nlohmann::json setup = _kernelSetups.front();
            _kernelSetups.pop();

            string sparsifierType = setup.value("sparsifier", "histogram_uniform");
            setup.erase("sparsifier");
            setup["sparse_param"] = setup.value("sparse_param", "density");

            for (const auto& [species, densities]: allDensitiesPerSpecies) {
                nlohmann::json setupPerSpecies = setup;

                if (setupPerSpecies.contains("species")) {
                    if (setupPerSpecies["species"].is_string()) {
                        setupPerSpecies["species"] = {setupPerSpecies["species"]};
                    }
                    if (!setupPerSpecies["species"].contains(species)) continue;
                    setupPerSpecies.erase("species");
                }

                auto sparsifier = ParserRegistry<Sparsifier>::getRegistry()[sparsifierType](setupPerSpecies);

                for (nlohmann::json& kernelParams : sparsifier->selectSparsePoints(densities)) {
                    kernelParams["species"] = species;
                    _kernels.push_back(ParserRegistry<EamKernel>::get(kernelParams));
                }
            }
        }

        mapKernelIds();
    }

    vector<Covariance> EamDescriptor::covariate(const AtomicStructure &atomicStructure) {
        vector<Covariance> result;

        EamKernelIndex kernelIndex = doIndex(atomicStructure);

        for (auto &[species, kernelIds]: _kernelIndicesPerSpecies) {
            if (!kernelIndex.contains(species)) kernelIndex[species] = {};

            for (auto& id: kernelIds) {
                result.push_back(_kernels[id]->covariance(atomicStructure, kernelIndex[species]));
            }
        }

        return result;
    }

    vector<shared_ptr<MatrixBlock>> EamDescriptor::selfCovariate() {

        vector<shared_ptr<MatrixBlock>> result;

        for (const auto &kernelIds: _kernelIndicesPerSpecies | views::values) {

            auto covariance = make_shared<MatrixBlock>(kernelIds.size(), kernelIds.size());

            for (size_t i = 0; i < kernelIds.size(); i++) {
                for (size_t j = i; j < kernelIds.size(); j++) {
                    (*covariance)(i, j) = _kernels[i]->crossCovariance(_kernels[j]);
                    (*covariance)(j, i) = (*covariance)(i, j);
                }
            }

            result.push_back(covariance);
        }

        return result;
    }

    Predictions EamDescriptor::predict(const AtomicStructure &atomicStructure) {
        auto indexes = doIndex(atomicStructure);
        Predictions result{};
        for (const auto& kernel: _kernels) {
            result = result + kernel->predict(atomicStructure, indexes[kernel->getFilter()]);
        }
        return result;
    }

    void EamDescriptor::tabulate(TabulationData &table) {
        auto& eamTable = table.newEamGrids();

        for (const auto& [species, kernelIds]: _kernelIndicesPerSpecies) {
            for (const auto& it: eamTable.getEnergyGrid(species)) {
                for (auto& id: kernelIds) {
                    it.value += _kernels[id]->value(it.pos) * _kernels[id]->coefficient.value();
                }
            }
        }

        for (const auto& contributorSpecies: _kernelIndicesPerSpecies | views::keys) {
            for (const auto& receiverSpecies: _kernelIndicesPerSpecies | views::keys) {

                auto speciesPair = ContributorReceiverSpecies{
                    .contributor = contributorSpecies, .receiver = receiverSpecies
                };

                auto pairFunction = _defaultPairFunction;
                if (_pairFunctions.contains(speciesPair)) {
                    pairFunction = _pairFunctions[speciesPair];
                }

                for (const auto& it: eamTable.getPairFunctionGrid(speciesPair)) {
                    it.value = pairFunction->evaluate(it.pos);
                }
            }
        }
    }

    EamKernelIndex EamDescriptor::doIndex(const AtomicStructure &structure) const {

        EamKernelIndex result{};

        for (size_t atomIdx = 0; atomIdx < structure.size(); atomIdx++) {

            double totalDensity = 0;
            vector<pair<NeighbourData, double>> densityDerivatives;

            Species species = structure.species[atomIdx];
            for (NeighbourData neighbour: structure.neighbours.value()[atomIdx]) {
                if (neighbour.distance > _maxCutoff) continue;

                ContributorReceiverSpecies orderedSpeciesPair = {
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

    void EamDescriptor::mapKernelIds() {
        _kernelIndicesPerSpecies.clear();
        for (size_t i = 0; i < _kernels.size(); i++) {
            if (!_kernelIndicesPerSpecies.contains(_kernels[i]->getFilter())) {
                _kernelIndicesPerSpecies[_kernels[i]->getFilter()] = {};
            }
            _kernelIndicesPerSpecies[_kernels[i]->getFilter()].push_back(i);
        }
    }
}
