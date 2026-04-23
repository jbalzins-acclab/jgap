#include "core/descriptors/eam/EamDescriptor.hpp"

#include <random>
#include <utility>

#include "../../kernels/eam/EamSE.hpp"
#include "core/descriptors/eam/pair_functions/FailOnUsePairFunction.hpp"
#include "io/log/StdoutLogger.hpp"
#include "io/parse/ParserRegistry.hpp"

namespace jgap {
    EamDescriptor::EamDescriptor(std::vector<std::shared_ptr<EamKernel>> kernels,
                                 std::shared_ptr<EamPairFunction> &default_pair_function,
                                 std::map<ContributorReceiverSpecies, std::shared_ptr<EamPairFunction>> pairFunctions)
        : kernels_(std::move(kernels)),
          pair_functions_(std::move(pairFunctions)) {

        if (default_pair_function == nullptr) {
            default_pair_function = std::make_shared<FailOnUsePairFunction>();
        } else {
            default_pair_function = default_pair_function_;
        }

        max_cutoff_ = 0;
        for (const auto& pf: pairFunctions | std::views::values) {
            max_cutoff_ = std::max(max_cutoff_, pf->getCutoff());
        }
        mapKernelIds();
    }

    std::shared_ptr<EamDescriptor> EamDescriptor::fromDataNode(const DataNode &params) {
        JGAP_LOG_DEBUG("Parsing EAM descriptor params");

        kernels_ = {};
        if (params.contains("kernels")) {
            for (const DataNode& kernelParams : params["kernels"].asArray()) {
                kernels_.push_back(REGISTRY_GET(EamKernel, kernelParams));
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

        max_cutoff_ = 0;
        for (const auto& pfParams: params["pair_functions"]) {

            auto pf = REGISTRY_GET(EamPairFunction, pfParams);

            if (pfParams.contains("species")) {
                auto s1 = pfParams["species"][0];
                auto s2 = pfParams["species"][1];
                pair_functions_[{s1, s2}] = pf;
                pair_functions_[{s2, s1}] = pf;
            } else if (pfParams.contains("contributor_receiver_species")) {
                auto s1 = pfParams["contributor_receiver_species"][0];
                auto s2 = pfParams["contributor_receiver_species"][1];
                pair_functions_[{s1, s2}] = pf;
            } else {
                default_pair_function_ = pf;
            }

            max_cutoff_ = std::max(max_cutoff_, pf->getCutoff());
        }
    }

    DataNode EamDescriptor::serialize() {

        auto kernelsData = DataNode::array();
        for (const auto &kernel : kernels_) {
            auto serial = std::dynamic_pointer_cast<Serializable>(kernel);

            kernelsData.pushBack(serial->serialize());
            kernelsData.back()["type"] = serial->getTypeId();
        }

        DataNode pfData = DataNode::array();
        auto defaultPf = std::dynamic_pointer_cast<Serializable>(default_pair_function_);

        auto defaultPairFunctionData = defaultPf->serialize();
        defaultPairFunctionData["type"] = defaultPf->getTypeId();
        pfData.pushBack(defaultPairFunctionData);

        for (const auto& [speciesPair, pf]: pair_functions_) {
            auto serial = std::dynamic_pointer_cast<Serializable>(pf);
            auto newPfData = serial->serialize();
            newPfData["type"] = serial->getTypeId();
            newPfData["contributor_receiver_species"] = std::vector{speciesPair.contributor, speciesPair.receiver};
            pfData.pushBack(newPfData);
        }

        return {
            {"kernels", kernelsData},
            {"pair_functions", pfData}
        };
    }

    CutoffRanges EamDescriptor::getCutoff() {
        auto res = CutoffRanges{.twoBody = max_cutoff_};

        double minDensity = 0.0, maxDensity = std::numeric_limits<double>::std::min();
        for (const auto& kernel: kernels_) {
            minDensity = std::min(minDensity, kernel->getDensityRange().first);
            maxDensity = std::max(maxDensity, kernel->getDensityRange().second);
        }
        res.minEam = minDensity;
        res.maxEam = maxDensity;

        return res;
    }

    std::vector<std::shared_ptr<IKernel>> EamDescriptor::getKernels() {
        std::vector<std::shared_ptr<IKernel>> res;
        for (const auto& kernelIds : kernel_indices_per_species_ | std::views::values) {
            for (const auto& id: kernelIds) {
                res.push_back(kernels_[id]);
            }
        }
        return res;
    }

    void EamDescriptor::setupSparseKernels(const std::vector<AtomicStructure> &fromData) {
        JGAP_LOG_INFO("Doing EAM sparsification from data");

        if (_kernelSetups.empty()) {
            JGAP_LOG_WARN("All 3b kernels were pre-set");
            return;
        }

        std::map<Species, std::vector<std::vector<double>>> allDensitiesPerSpecies;
        std::vector<EamKernelIndex> indexArr;
        for (const auto& structure : fromData) {
            for (auto structureIndex = doIndex(structure);
                 const auto& [species, densities]: structureIndex) {
                if (!allDensitiesPerSpecies.contains(species)) {
                    allDensitiesPerSpecies[species] = {};
                }
                for (const auto& densityData: densities) {
                    // TODO: low density cap - avoid zero of isolated_atom density(?):
                    // if (densityData.density < 0.2) continue;
                    allDensitiesPerSpecies[species].push_back(std::vector{densityData.density});
                }
            }
        }

        while (!_kernelSetups.empty()) {
            DataNode setup = _kernelSetups.front();
            _kernelSetups.pop();

            std::string sparsifierType = setup.value("sparsifier", "histogram_uniform");
            setup.erase("sparsifier");
            setup["sparse_param"] = setup.value("sparse_param", "density");

            for (const auto& [species, densities]: allDensitiesPerSpecies) {
                DataNode setupPerSpecies = setup;

                if (setupPerSpecies.contains("species")) {
                    if (setupPerSpecies["species"].is_string()) {
                        setupPerSpecies["species"] = {setupPerSpecies["species"]};
                    }
                    if (!setupPerSpecies["species"].contains(species)) continue;
                    setupPerSpecies.erase("species");
                }

                auto sparsifier = ParserRegistry<Sparsifier>::getRegistry()[sparsifierType](setupPerSpecies);

                for (DataNode& kernelParams : sparsifier->selectSparsePoints(densities)) {
                    kernelParams["species"] = species;
                    kernels_.push_back(ParserRegistry<EamKernel>::get(kernelParams));
                }
            }
        }

        mapKernelIds();
    }

    std::vector<Covariance> EamDescriptor::covariate(const AtomicStructure &atomicStructure) {
        std::vector<Covariance> result;

        EamKernelIndex kernelIndex = doIndex(atomicStructure);

        for (auto &[species, kernelIds]: kernel_indices_per_species_) {
            if (!kernelIndex.contains(species)) kernelIndex[species] = {};

            for (auto& id: kernelIds) {
                result.push_back(kernels_[id]->covariance(atomicStructure, kernelIndex[species]));
            }
        }

        return result;
    }

    std::vector<std::shared_ptr<MatrixBlock>> EamDescriptor::selfCovariate() {

        std::vector<std::shared_ptr<MatrixBlock>> result;

        for (const auto &kernelIds: kernel_indices_per_species_ | std::views::values) {

            auto covariance = std::make_shared<MatrixBlock>(kernelIds.size(), kernelIds.size());

            for (size_t i = 0; i < kernelIds.size(); i++) {
                for (size_t j = i; j < kernelIds.size(); j++) {
                    (*covariance)(i, j) = kernels_[i]->crossCovariance(kernels_[j]);
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
        for (const auto& kernel: kernels_) {
            result = result + kernel->predict(atomicStructure, indexes[kernel->getFilter()]);
        }
        return result;
    }

    void EamDescriptor::tabulate(TabulationData &table) {
        auto& eamTable = table.newEamGrids();

        for (const auto& [species, kernelIds]: kernel_indices_per_species_) {
            for (const auto& it: eamTable.getOrMakeEnergyGrid(species)) {
                for (auto& id: kernelIds) {
                    it.value += kernels_[id]->value(it.pos) * kernels_[id]->coefficient.value();
                }
            }
        }

        for (const auto& contributorSpecies: kernel_indices_per_species_ | std::views::keys) {
            for (const auto& receiverSpecies: kernel_indices_per_species_ | std::views::keys) {

                auto speciesPair = ContributorReceiverSpecies{
                    .contributor = contributorSpecies, .receiver = receiverSpecies
                };

                auto pairFunction = default_pair_function_;
                if (pair_functions_.contains(speciesPair)) {
                    pairFunction = pair_functions_[speciesPair];
                }

                for (const auto& it: eamTable.getOrMakePairFunctionGrid(speciesPair)) {
                    it.value = pairFunction->evaluate(it.pos);
                }
            }
        }
    }

    EamKernelIndex EamDescriptor::doIndex(const AtomicStructure &structure) const {

        EamKernelIndex result{};

        for (size_t atomIdx = 0; atomIdx < structure.size(); atomIdx++) {

            double totalDensity = 0;
            std::vector<std::pair<NeighbourData, double>> densityDerivatives;

            Species species = structure.species[atomIdx];
            for (NeighbourData neighbour: structure.neighbours_ascending_separation.value()[atomIdx]) {
                if (neighbour.distance > max_cutoff_) continue;

                ContributorReceiverSpecies orderedSpeciesPair = {
                    structure.species[neighbour.index],  species
                };

                std::shared_ptr<EamPairFunction> pf = default_pair_function_;
                if (pair_functions_.contains(orderedSpeciesPair)) {
                    pf = pair_functions_.at(orderedSpeciesPair);
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
        kernel_indices_per_species_.clear();
        for (size_t i = 0; i < kernels_.size(); i++) {
            if (!kernel_indices_per_species_.contains(kernels_[i]->getFilter())) {
                kernel_indices_per_species_[kernels_[i]->getFilter()] = {};
            }
            kernel_indices_per_species_[kernels_[i]->getFilter()].push_back(i);
        }
    }
}
