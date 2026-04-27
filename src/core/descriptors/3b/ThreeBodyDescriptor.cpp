#include "core/descriptors/3b/ThreeBodyDescriptor.hpp"

#include <random>
#include <tbb/parallel_for_each.h>

#include "core/descriptors/3b/ThreeBodySE.hpp"
#include "io/log/StdoutLogger.hpp"
#include "utils/Utils.hpp"

namespace jgap {
    ThreeBodyDescriptor::ThreeBodyDescriptor(shared_ptr<CutoffFunction> &cutoffFunction,
                                             vector<shared_ptr<ThreeBodyKernel>> &kernels)
        : _cutoff(cutoffFunction->getCutoff()),
          _cutoffFunction(std::move(cutoffFunction)),
          _kernels(std::move(kernels)) {

        mapKernelIds();
    }

    ThreeBodyDescriptor::ThreeBodyDescriptor(const nlohmann::json &params) {
        CurrentLogger::get()->debug("Parsing 3b descriptor params");

        _cutoffFunction = ParserRegistry<CutoffFunction>::get(require(params, "cutoff"));
        _cutoff = _cutoffFunction->getCutoff();

        _kernels = {};
        if (params.contains("kernels")) {
            for (const nlohmann::json& kernelParams : params["kernels"]) {
                _kernels.push_back(ParserRegistry<ThreeBodyKernel>::get(kernelParams));
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
    }

    nlohmann::json ThreeBodyDescriptor::serialize() {
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

    void ThreeBodyDescriptor::setupKernels(const vector<AtomicStructure> &fromData) {
        CurrentLogger::get()->info("Doing 3b sparsification from data");

        if (_kernelSetups.empty()) {
            CurrentLogger::get()->warn("All 3b kernels were pre-set");
            return;
        }

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

        while (!_kernelSetups.empty()) {
            nlohmann::json setup = _kernelSetups.front();
            _kernelSetups.pop();

            string sparsifierType = setup.value("sparsifier", "histogram_uniform");
            setup.erase("sparsifier");
            setup["sparse_param"] = setup.value("sparse_param", "q");

            for (const auto &[tripletInData, tripletsPerSpecies]: allTriplets) {
                nlohmann::json setupPerSpecies = setup;

                if (!checkSpecies(tripletInData, setupPerSpecies.value("species", nlohmann::json::array()))) {
                    continue;
                }

                setupPerSpecies.erase("species");
                auto sparsifier = ParserRegistry<Sparsifier>::getRegistry()[sparsifierType](setupPerSpecies);

                for (nlohmann::json& kernelParams : sparsifier->selectSparsePoints(tripletsPerSpecies)) {
                    kernelParams["species_triplet"] = vector{
                        tripletInData.root, tripletInData.nodes.first(), tripletInData.nodes.second()
                    };
                    kernelParams["descriptor_prefactors"] = invariantTripletToCutoff({
                        kernelParams["q"][0], kernelParams["q"][1], kernelParams["q"][2]
                    });
                    _kernels.push_back(ParserRegistry<ThreeBodyKernel>::get(kernelParams));
                }
            }
        }
        mapKernelIds();
    }

    vector<shared_ptr<IKernel>> ThreeBodyDescriptor::getKernels() {
        vector<shared_ptr<IKernel>> res;
        for (auto& kernelIds : _kernelIdsPerSpeciesTriplet | views::values) {
            for (const auto& kernelId : kernelIds) {
                res.push_back(static_pointer_cast<IKernel>(_kernels[kernelId]));
            }
        }
        return res;
    }

    PotentialPrediction ThreeBodyDescriptor::predict(const AtomicStructure &atomicStructure) {
        auto indexes = doIndex(atomicStructure);
        PotentialPrediction result{};
        for (const auto& kernel: _kernels) {
            result = result + kernel->predict(atomicStructure, indexes[kernel->getFilter()]);
        }
        return result;
    }

    vector<Covariance> ThreeBodyDescriptor::covariate(const AtomicStructure &atomicStructure) {

        auto indexes = doIndex(atomicStructure);

        auto covariates = vector<Covariance>();
        for (auto &kernelIndices: _kernelIdsPerSpeciesTriplet | views::values) {
            for (const auto &kernelId: kernelIndices) {
                auto& kernel = _kernels[kernelId];
                covariates.push_back(
                    kernel->covariance(atomicStructure, GET_OR_DEFAULT(indexes, kernel->getFilter(), ThreeBodyIndex{}))
                );
            }
        }

        return covariates;
    }

    vector<shared_ptr<MatrixBlock>> ThreeBodyDescriptor::selfCovariate() {
        vector<shared_ptr<MatrixBlock>> result;

        for (auto &kernelIds: _kernelIdsPerSpeciesTriplet | views::values) {

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

    TabulationData ThreeBodyDescriptor::tabulate(const TabulationParams &params) {

        TabulationData result{};

        for (const auto &[speciesTriplet, kernelIds]: _kernelIdsPerSpeciesTriplet) {

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

                for (const size_t kernelId: kernelIds) {
                    tripletEnergies[iGrid[0]][iGrid[1]][iGrid[2]] += _kernels[kernelId]->value(
                        {.q = invariantTriplet, .fCut = invariantTripletToCutoff(invariantTriplet)}
                        ) * 2.0/*q_ijk + q_jik*/ * _kernels[kernelId]->coefficient.value();
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

    Vector3 ThreeBodyDescriptor::toInvariantTriplet(const double r01, const double r02, const double r12) {
        return {r01 + r02, (r01-r02) * (r01 - r02), r12};
    }

    array<Vector3, 3> ThreeBodyDescriptor::invariantTripletGradients(const double r01, const double r02) {
        return {
            Vector3{1, 2 * (r01 - r02), 0},
            Vector3{1, 2 * (r02 - r01), 0},
            Vector3{0, 0, 1},
        };
    }

    map<SpeciesTriplet, ThreeBodyIndex> ThreeBodyDescriptor::doIndex(
                                                const AtomicStructure &atomicStructure) const {

        map<SpeciesTriplet, ThreeBodyIndex> indexes;

        for (size_t atomIndex = 0; atomIndex < atomicStructure.size(); atomIndex++) {
            auto atom0 = atomicStructure[atomIndex];

            // "nl" = neighbourList

            vector<size_t> usefulNLIndexes{};
            usefulNLIndexes.reserve(atom0.neighbours().size());

            for (size_t nlIndex1 = 0; nlIndex1 < atom0.neighbours().size(); nlIndex1++) {
                auto neighbour1 = atom0.neighbours()[nlIndex1];
                auto atom1 = atomicStructure[neighbour1.index];
                if (neighbour1.distance > _cutoff) continue;

                for (const size_t& nlIndex2: usefulNLIndexes) {
                    auto neighbour2 = atom0.neighbours()[nlIndex2];
                    auto atom2 = atomicStructure[neighbour2.index];

                    auto speciesTriplet = SpeciesTriplet{atom0.species(),{atom1.species(), atom2.species()}};

                    if (!indexes.contains(speciesTriplet)) {
                        indexes[speciesTriplet] = ThreeBodyIndex();
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

                // !!! ORDER SENSITIVE !!!
                usefulNLIndexes.emplace_back(nlIndex1);
            }
        }

        return indexes;
    }

    bool ThreeBodyDescriptor::checkSpecies(const SpeciesTriplet& tripletInData, nlohmann::json filters) {

        if (!filters.is_array()) {
            CurrentLogger::get()->logAndThrow("3b species filter is non-array: {}", filters.dump());
        }
        if (filters.empty()) return true;

        if (filters[0].is_string()) {
            // "species": ["Fe", "Ni", "Cr"] => "species": [["Fe", "Ni", "Cr"]]
            filters = nlohmann::json::array({filters});
        }

        bool passedAFilter = false;
        for (const auto &filter: filters) {
            if (!filter.is_array() || filter.size() != 3) {
                CurrentLogger::get()->logAndThrow("Wrong 3b species specs: {}", filters.dump());
            }
            if (tripletInData == SpeciesTriplet(filter[0], {filter[1], filter[2]})) {
                passedAFilter = true;
                break;
            }
        }
        return passedAFilter;
    }

    void ThreeBodyDescriptor::mapKernelIds() {
        _kernelIdsPerSpeciesTriplet.clear();
        for (size_t i = 0; i < _kernels.size(); i++) {
            if (!_kernelIdsPerSpeciesTriplet.contains(_kernels[i]->getFilter())) {
                _kernelIdsPerSpeciesTriplet[_kernels[i]->getFilter()] = {};
            }
            _kernelIdsPerSpeciesTriplet[_kernels[i]->getFilter()].push_back(i);
        }
    }
}
