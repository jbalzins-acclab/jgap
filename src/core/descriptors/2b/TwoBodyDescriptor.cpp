#include "core/descriptors/2b/TwoBodyDescriptor.hpp"

#include "../../kernels/2b/TwoBodySE.hpp"
#include "io/log/StdoutLogger.hpp"
#include "io/parse/ParserRegistry.hpp"
#include "utils/Utils.hpp"

namespace jgap {
    TwoBodyDescriptor::TwoBodyDescriptor(std::shared_ptr<CutoffFunction> cutoffFunction,
                                         std::vector<std::shared_ptr<TwoBodyKernel>> kernels)
        : _cutoffFunction(std::move(cutoffFunction)), _kernels(std::move(kernels)) {
        mapKernelIds();
    }

    TwoBodyDescriptor::TwoBodyDescriptor(const DataNode& params) {
        JGAP_LOG_DEBUG("Parsing 2b descriptor params");

        _cutoffFunction = ParserRegistry<CutoffFunction>::get(params.["cutoff"]);

        _kernels = {};
        if (params.contains("kernels")) {
            for (const DataNode& kernelParams : params["kernels"]) {
                _kernels.push_back(ParserRegistry<TwoBodyKernel>::get(kernelParams));
            }
        }
        mapKernelIds();

        _kernelSetups = {};
        if (params.contains("kernel_setups")) {
            for (const DataNode& kernelSetup : params["kernel_setups"]) {
                _kernelSetups.push(kernelSetup);
            }
        }
        if (params.contains("kernel_setup")) {
            _kernelSetups.push(params["kernel_setup"]);
        }

        if (_kernelSetups.empty() && _kernels.empty()) {
            JGAP_LOG_AND_THROW("Descriptor with no kernels");
        }
    }

    DataNode TwoBodyDescriptor::serialize() {

        auto kernelsData = DataNode::array();
        for (const auto &kernel : _kernels) {
            kernelsData.pushBack(kernel->serialize());
            kernelsData.back()["type"] = kernel->getType();
        }

        auto cutoffData = _cutoffFunction->serialize();
        cutoffData["type"] = _cutoffFunction->getType();

        return {
            {"kernels", kernelsData},
            {"cutoff", cutoffData}
        };
    }

    std::vector<std::shared_ptr<IKernel>> TwoBodyDescriptor::getKernels() {
        std::vector<std::shared_ptr<IKernel>> res;
        for (auto& kernelIds : _kernelIdsPerSpeciesPair | std::views::values) {
            for (const auto& kernelId : kernelIds) {
                res.push_back(static_pointer_cast<IKernel>(_kernels[kernelId]));
            }
        }
        return res;
    }

    void TwoBodyDescriptor::setupSparseKernels(const std::vector<AtomicStructure> &fromData) {
        JGAP_LOG_INFO("Doing 2b sparsification from data");

        if (_kernelSetups.empty()) {
            JGAP_LOG_WARN("All 2b kernels were pre-std::set");
            return;
        }

        std::map<SpeciesPair, std::vector<std::vector<double>>> allDistances;
        for (const auto &structure: fromData) {
            const auto structureIndex = doIndex(structure);

            for (const auto &[speciesPair, perPairIndex]: structureIndex) {
                if (!allDistances.contains(speciesPair)) allDistances[speciesPair] = {};
                for (const auto &entity: perPairIndex) {
                    allDistances[speciesPair].push_back(std::vector{entity.r});
                }
            }
        }

        while (!_kernelSetups.empty()) {
            DataNode setup = _kernelSetups.front();
            _kernelSetups.pop();

            std::string sparsifierType = setup.value("sparsifier", "histogram_uniform");
            setup.erase("sparsifier");
            setup["sparse_param"] = setup.value("sparse_param", "r");

            for (const auto &[pairInData, distancesPerPair]: allDistances) {
                DataNode setupPerPair = setup; // modified for the specific std::pair

                if (!checkSpecies(pairInData, setupPerPair.value("species", DataNode::array()))) {
                    continue;
                }

                setupPerPair.erase("species");
                auto sparsifier = ParserRegistry<Sparsifier>::getRegistry()[sparsifierType](setupPerPair);

                for (DataNode& kernelParams : sparsifier->selectSparsePoints(distancesPerPair)) {
                    kernelParams["species_pair"] = std::vector{pairInData.first(), pairInData.second()};
                    kernelParams["descriptor_prefactors"] = _cutoffFunction->evaluate(kernelParams["r"]);
                    _kernels.push_back(ParserRegistry<TwoBodyKernel>::get(kernelParams));
                }
            }
        }
        mapKernelIds();
    }

    std::vector<Covariance> TwoBodyDescriptor::covariate(const AtomicStructure &atomicStructure) {

        auto indexes = doIndex(atomicStructure);

        auto covariates = std::vector<Covariance>();
        for (const auto &kernel: _kernels) {
            covariates.push_back(
                kernel->covariance(atomicStructure, GET_OR_DEFAULT(indexes, kernel->getFilter(), TwoBodyIndex{}))
                );
        }

        return covariates;
    }

    std::vector<std::shared_ptr<MatrixBlock>> TwoBodyDescriptor::selfCovariate() {
        std::vector<std::shared_ptr<MatrixBlock>> result;

        for (auto &kernelIndices: _kernelIdsPerSpeciesPair | std::views::values) {

            auto covariance = std::make_shared<MatrixBlock>(kernelIndices.size(), kernelIndices.size());

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

    Predictions TwoBodyDescriptor::predict(const AtomicStructure &atomicStructure) {
        auto indexes = doIndex(atomicStructure);
        Predictions result{};
        for (const auto& kernel: _kernels) {
            result = result + kernel->predict(atomicStructure, indexes[kernel->getFilter()]);
        }
        return result;
    }

    void TwoBodyDescriptor::tabulate(TabulationData &table) {

        for (const auto &[speciesPair, kernelIds]: _kernelIdsPerSpeciesPair) {
            for (const auto& it: table.getOrMake2bGrid(speciesPair)) {
                for (const size_t kernelId: kernelIds) {
                    it.value += _kernels[kernelId]->value( {
                        .r=it.pos, .fCut = _cutoffFunction->evaluate(it.pos)
                        }) * _kernels[kernelId]->coefficient.value() * 2.0/*K_ij+K_ji*/;
                }
            }
        }
    }

    std::map<SpeciesPair, TwoBodyIndex> TwoBodyDescriptor::doIndex(const AtomicStructure &atomicStructure) const {

        std::map<SpeciesPair, TwoBodyIndex> indexes;

        for (size_t atomIndex = 0; atomIndex < atomicStructure.size(); atomIndex++) {
            auto atom = atomicStructure[atomIndex];

            for (size_t neighbourListIndex = 0; neighbourListIndex < atom.neighboursAscendingSeparation().size(); neighbourListIndex++) {
                const auto neighbour = atom.neighboursAscendingSeparation()[neighbourListIndex];

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

    bool TwoBodyDescriptor::checkSpecies(const SpeciesPair& pairInData, DataNode filters) {
        if (!filters()) {
            JGAP_LOG_AND_THROW("2b species filter is non-array: {}", filters.dump());
        }

        if (filters.empty()) return true;

        if (filters[0].is_string()) {
            // "species": ["Fe", "Ni"] => "species": [["Fe", "Ni"]]
            filters = DataNode({filters});
        }

        bool passedAFilter = false;
        for (const auto &filter: filters) {
            if (!filter.is_array() || filter.size() != 2) {
                JGAP_LOG_AND_THROW("Wrong 2b species specs: {}", filters.dump());
            }
            if (pairInData == SpeciesPair{filter[0], filter[1]}) {
                passedAFilter = true;
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
