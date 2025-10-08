#include "core/descriptors/2b/TwoBodyDescriptor.hpp"

#include <nlohmann/json.hpp>

#include "core/descriptors/2b/TwoBodySE.hpp"
#include "io/log/StdoutLogger.hpp"
#include "io/parse/ParserRegistry.hpp"
#include "utils/Utils.hpp"

using namespace std;

namespace jgap {
    TwoBodyDescriptor::TwoBodyDescriptor(shared_ptr<CutoffFunction> cutoffFunction,
                                         vector<shared_ptr<TwoBodyKernel>> kernels)
        : _cutoffFunction(std::move(cutoffFunction)), _kernels(std::move(kernels)) {
        mapKernelIds();
    }

    TwoBodyDescriptor::TwoBodyDescriptor(const nlohmann::json& params) {
        CurrentLogger::get()->debug("Parsing 2b descriptor params");

        _cutoffFunction = ParserRegistry<CutoffFunction>::get(params["cutoff"]);

        _kernels = {};
        if (params.contains("kernels")) {
            for (const nlohmann::json& kernelParams : params["kernels"]) {
                _kernels.push_back(ParserRegistry<TwoBodyKernel>::get(kernelParams));
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

        if (_kernelSetups.empty() && _kernels.empty()) {
            CurrentLogger::get()->logAndThrow("Descriptor with no kernels");
        }
    }

    nlohmann::json TwoBodyDescriptor::serialize() {

        auto kernelsData = nlohmann::json::array();
        for (const auto &kernel : _kernels) {
            kernelsData.push_back(kernel->serialize());
            kernelsData.back()["type"] = kernel->getType();
        }

        auto cutoffData = _cutoffFunction->serialize();
        cutoffData["type"] = _cutoffFunction->getType();

        return {
            {"kernels", kernelsData},
            {"cutoff", cutoffData}
        };
    }

    vector<shared_ptr<IKernel>> TwoBodyDescriptor::getKernels() {
        vector<shared_ptr<IKernel>> res;
        for (auto& kernelIds : _kernelIdsPerSpeciesPair | views::values) {
            for (const auto& kernelId : kernelIds) {
                res.push_back(static_pointer_cast<IKernel>(_kernels[kernelId]));
            }
        }
        return res;
    }

    void TwoBodyDescriptor::setupKernels(const vector<AtomicStructure> &fromData) {
        CurrentLogger::get()->info("Doing 2b sparsification from data");

        if (_kernelSetups.empty()) {
            CurrentLogger::get()->warn("All 2b kenels were pre-set");
            return;
        }

        map<SpeciesPair, vector<vector<double>>> allDistances;
        for (const auto &structure: fromData) {
            const auto structureIndex = doIndex(structure);

            for (const auto &[speciesPair, perPairIndex]: structureIndex) {
                if (!allDistances.contains(speciesPair)) allDistances[speciesPair] = {};
                for (const auto &entity: perPairIndex) {
                    allDistances[speciesPair].push_back(vector{entity.r});
                }
            }
        }

        while (!_kernelSetups.empty()) {
            nlohmann::json setup = _kernelSetups.front();
            _kernelSetups.pop();

            string sparsifierType = setup.value("sparsifier", "histogram_uniform");
            setup.erase("sparsifier");
            setup["sparse_param"] = setup.value("sparse_param", "r");

            for (const auto &[pairInData, distancesPerPair]: allDistances) {
                nlohmann::json setupPerPair = setup; // modified for the specific pair

                if (!checkSpecies(pairInData, setupPerPair.value("species", nlohmann::json::array()))) {
                    continue;
                }

                setupPerPair.erase("species");
                auto sparsifier = ParserRegistry<Sparsifier>::getRegistry()[sparsifierType](setupPerPair);

                for (nlohmann::json& kernelParams : sparsifier->selectSparsePoints(distancesPerPair)) {
                    kernelParams["species_pair"] = vector{pairInData.first(), pairInData.second()};
                    kernelParams["descriptor_prefactors"] = _cutoffFunction->evaluate(kernelParams["r"]);
                    _kernels.push_back(ParserRegistry<TwoBodyKernel>::get(kernelParams));
                }
            }
        }
        mapKernelIds();
    }

    vector<Covariance> TwoBodyDescriptor::covariate(const AtomicStructure &atomicStructure) {

        auto indexes = doIndex(atomicStructure);

        auto covariates = vector<Covariance>();
        for (const auto &kernel: _kernels) {
            covariates.push_back(
                kernel->covariance(atomicStructure, GET_OR_DEFAULT(indexes, kernel->getFilter(), TwoBodyIndex{}))
                );
        }

        return covariates;
    }

    vector<shared_ptr<MatrixBlock>> TwoBodyDescriptor::selfCovariate() {
        vector<shared_ptr<MatrixBlock>> result;

        for (auto &kernelIndices: _kernelIdsPerSpeciesPair | views::values) {

            auto covariance = make_shared<MatrixBlock>(kernelIndices.size(), kernelIndices.size());

            for (size_t i = 0; i < kernelIndices.size(); i++) {
                for (size_t j = i; j < kernelIndices.size(); j++) {
                    (*covariance)(i, j) = _kernels[i]->crossCovariance(_kernels[j]);
                    (*covariance)(j, i) = (*covariance)(i, j);
                }
            }

            result.push_back(covariance);
        }

        return result;
    }

    PotentialPrediction TwoBodyDescriptor::predict(const AtomicStructure &atomicStructure) {
        auto indexes = doIndex(atomicStructure);
        PotentialPrediction result{};
        for (const auto& kernel: _kernels) {
            result = result + kernel->predict(atomicStructure, indexes[kernel->getFilter()]);
        }
        return result;
    }

    TabulationData TwoBodyDescriptor::tabulate(const TabulationParams &params) {

        TabulationData result{};

        for (const auto &[speciesPair, kernelIds]: _kernelIdsPerSpeciesPair) {

            auto pairEnergies = vector(params.grid2b.size(), 0.0);

            for (size_t iGrid = 0; iGrid < params.grid2b.size(); iGrid++) {
                for (size_t kernelId: kernelIds) {
                    pairEnergies[iGrid] += _kernels[kernelId]->value( {
                        .r=params.grid2b[iGrid], .fCut = _cutoffFunction->evaluate(params.grid2b[iGrid])
                    }) * _kernels[kernelId]->coefficient.value() * 2.0/*K_ij+K_ji*/;
                }
            }

            result.pairEnergies[speciesPair] = pairEnergies;
        }

        return result;
    }

    map<SpeciesPair, TwoBodyIndex> TwoBodyDescriptor::doIndex(const AtomicStructure &atomicStructure) const {

        map<SpeciesPair, TwoBodyIndex> indexes;

        for (size_t atomIndex = 0; atomIndex < atomicStructure.size(); atomIndex++) {
            auto atom = atomicStructure[atomIndex];

            for (size_t neighbourListIndex = 0; neighbourListIndex < atom.neighbours().size(); neighbourListIndex++) {
                const auto neighbour = atom.neighbours()[neighbourListIndex];

                if (neighbour.index < atomIndex) continue;
                if (neighbour.distance > _cutoffFunction->getCutoff()) continue;

                auto speciesPair = SpeciesPair{atom.species(), atomicStructure.species[neighbour.index]};
                if (!indexes.contains(speciesPair)) {
                    indexes[speciesPair] = TwoBodyIndex();
                }

                indexes[speciesPair].push_back({
                    .atomIndex0 = atomIndex,
                    .atomIndex1 = neighbour.index,
                    .r01 = atomicStructure.positions[neighbour.index] + neighbour.offset - atom.position(),
                    .r = neighbour.distance,
                    .fCut = _cutoffFunction->evaluate(neighbour.distance),
                    .dCut_dr = _cutoffFunction->differentiate(neighbour.distance)
                });
            }
        }

        return indexes;
    }

    bool TwoBodyDescriptor::checkSpecies(SpeciesPair pairInData, nlohmann::json filters) {
        if (!filters.is_array()) {
            CurrentLogger::get()->logAndThrow("2b species filter is non-array: {}", filters.dump());
        }
        CurrentLogger::get()->debug("{}:{}", pairInData.toString(), filters.dump());
        if (filters.empty()) return true;
        CurrentLogger::get()->debug("xx{}:{}", pairInData.toString(), filters.dump());

        if (filters[0].is_string()) {
            // "species": ["Fe", "Ni"] => "species": [["Fe", "Ni"]]
            filters = nlohmann::json::array({filters});
        }

        bool passedAFilter = false;
        for (const auto &filter: filters) {
            if (!filter.is_array() || filter.size() != 2) {
                CurrentLogger::get()->logAndThrow("Wrong 2b species specs: {}", filters.dump());
            }
            if (pairInData == SpeciesPair{filter[0], filter[1]}) {
                passedAFilter = true;
                CurrentLogger::get()->debug("xxxx{}:{}", pairInData.toString(), filters.dump());
                break;
            }
        }
        return passedAFilter;
    }

    void TwoBodyDescriptor::mapKernelIds() {
        _kernelIdsPerSpeciesPair.clear();
        for (size_t i = 0; i < _kernels.size(); i++) {
            if (!_kernelIdsPerSpeciesPair.contains(_kernels[i]->getFilter())) {
                _kernelIdsPerSpeciesPair[_kernels[i]->getFilter()] = {};
            }
            _kernelIdsPerSpeciesPair[_kernels[i]->getFilter()].push_back(i);
        }
    }
}
