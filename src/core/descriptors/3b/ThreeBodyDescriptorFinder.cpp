#include "core/descriptors/3b/ThreeBodyDescriptorFinder.hpp"

#include <random>
#include <tbb/parallel_for_each.h>

#include "io/log/StdoutLogger.hpp"
#include "utils/Utils.hpp"

namespace jgap {

    DescriptorFinder<3, 3>::DescriptorsFiltered ThreeBodyDescriptorFinder::
    find(const AtomicStructure &atomic_structure) {

    }

    DataNode ThreeBodyDescriptorFinder::serialize() {
        auto kernelsData = DataNode::array();
        for (const auto &kernel : _kernels) {
            kernelsData.pushBack(kernel->serialize());
            kernelsData.back()["type"] = kernel->getType();
        }

        auto cutoffData = cutoff_function_->serialize();
        cutoffData["type"] = cutoff_function_->getType();

        return {
                {"kernels", kernelsData},
                {"cutoff", cutoffData}
        };
    }

    void ThreeBodyDescriptorFinder::setupSparseKernels(const std::vector<AtomicStructure> &fromData) {
        JGAP_LOG_INFO("Doing 3b sparsification from data");

        if (_kernelSetups.empty()) {
            JGAP_LOG_WARN("All 3b kernels were pre-std::set");
            return;
        }

        std::map<SpeciesTriplet, std::vector<std::vector<double>>> allTriplets;
        for (const auto &structure: fromData) {
            for (const auto structureTriplets = doIndex(structure);
                 const auto &[speciesTriplet, points]: structureTriplets) {
                if (!allTriplets.contains(speciesTriplet)) {
                    allTriplets[speciesTriplet] = {};
                }

                for (const auto &point: points) {
                    allTriplets[speciesTriplet].push_back(std::vector{point.q.x, point.q.y, point.q.z});
                }
            }
        }

        while (!_kernelSetups.empty()) {
            DataNode setup = _kernelSetups.front();
            _kernelSetups.pop();

            std::string sparsifierType = setup.value("sparsifier", "histogram_uniform");
            setup.erase("sparsifier");
            setup["sparse_param"] = setup.value("sparse_param", "q");

            for (const auto &[tripletInData, tripletsPerSpecies]: allTriplets) {
                DataNode setupPerSpecies = setup;

                if (!checkSpecies(tripletInData, setupPerSpecies.value("species", DataNode::array()))) {
                    continue;
                }

                setupPerSpecies.erase("species");
                auto sparsifier = ParserRegistry<Sparsifier>::getRegistry()[sparsifierType](setupPerSpecies);

                for (DataNode& kernelParams : sparsifier->selectSparsePoints(tripletsPerSpecies)) {
                    kernelParams["species_triplet"] = std::vector{
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

    Predictions ThreeBodyDescriptorFinder::predict(const AtomicStructure &atomicStructure) {
        auto indexes = doIndex(atomicStructure);
        Predictions result{};
        for (const auto& kernel: _kernels) {
            result = result + kernel->predict(atomicStructure, indexes[kernel->getFilter()]);
        }
        return result;
    }

    std::vector<Covariance> ThreeBodyDescriptorFinder::covariate(const AtomicStructure &atomicStructure) {

        auto indexes = doIndex(atomicStructure);

        auto covariates = std::vector<Covariance>();
        for (const auto &kernel: _kernels) {
            covariates.push_back(
                kernel->covariance(atomicStructure, GET_OR_DEFAULT(indexes, kernel->getFilter(), ThreeBodyIndex{}))
                );
        }

        return covariates;
    }

    std::vector<std::shared_ptr<MatrixBlock>> ThreeBodyDescriptorFinder::selfCovariate() {
        std::vector<std::shared_ptr<MatrixBlock>> result;

        for (auto &kernelIds: _kernelIdsPerSpeciesTriplet | std::views::values) {

            auto covariance = std::make_shared<MatrixBlock>(kernelIds.size(), kernelIds.size());

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

    void ThreeBodyDescriptorFinder::tabulate(TabulationData &table) {

        for (const auto &[speciesTriplet, kernelIds]: _kernelIdsPerSpeciesTriplet) {

            auto& grid = table.getOrMake3bGrid(speciesTriplet);

            tbb::parallel_for_each(grid.begin(), grid.end(), [&](Grid3d::CellRef&& iter) {

                if (iter.index[0] > iter.index[1]) return;

                const Vector3 invariantTriplet = toInvariantTriplet(
                    iter.pos.x,
                    iter.pos.y,
                    sqrt(std::max/*numeric safety*/(
                        pow(iter.pos.x, 2) + pow(iter.pos.y, 2) - 2.0 * iter.pos.x * iter.pos.y * iter.pos.z, 0.0
                        ))
                );

                double contribution = 0.0;
                for (const size_t kernelId: kernelIds) {
                    contribution += _kernels[kernelId]->value(
                        {.q = invariantTriplet, .f_cut = invariantTripletToCutoff(invariantTriplet)}
                        ) * 2.0/*q_ijk + q_jik*/ * _kernels[kernelId]->coefficient.value();
                }
                iter.value += contribution;
                if (iter.index[0] != iter.index[1]) {
                    grid(iter.index[1], iter.index[0], iter.index[2]) += contribution;
                }
            });
        }

    }

    double ThreeBodyDescriptorFinder::invariantTripletToCutoff(const Vector3 &q) const {
        const double dDiff = sqrt(q.y);
        const double d1 = (dDiff + q.x) / 2.0;
        const double d2 = (q.x - dDiff) / 2.0;

        return cutoff_function_->evaluate(d1) * cutoff_function_->evaluate(d2);
    }

    Vector3 ThreeBodyDescriptorFinder::toInvariantTriplet(const double r01, const double r02, const double r12) const {
        return {r01 + r02, (r01-r02) * (r01 - r02), r12};
    }

    std::array<Vector3, 3> ThreeBodyDescriptorFinder::invariantTripletGradients(const double r01, const double r02) const {
        return {
            Vector3{1, 2 * (r01 - r02), 0},
            Vector3{1, 2 * (r02 - r01), 0},
            Vector3{0, 0, 1},
        };
    }

    std::vector<std::shared_ptr<IKernel>> ThreeBodyKernelCollection::getKernels() {
    }

    std::map<SpeciesTriplet, ThreeBodyIndex> ThreeBodyDescriptorFinder::doIndex(
                                                const AtomicStructure &atomic_structure) const {

        std::map<SpeciesTriplet, ThreeBodyIndex> indexes;

        const double cutoff = cutoff_function_->getCutoff();

        for (size_t atomIndex = 0; atomIndex < atomic_structure.size(); atomIndex++) {
            auto atom0 = atomic_structure[atomIndex];

            for (size_t nlIndex1 = 0; nlIndex1 < atom0.neighboursAscendingSeparation().size(); nlIndex1++) {
                auto neighbour1 = atom0.neighboursAscendingSeparation()[nlIndex1];
                auto atom1 = atomic_structure[neighbour1.index];
                if (neighbour1.distance > cutoff) continue;

                for (const auto& neighbour2: atom0.neighboursAscendingSeparation()) {
                    auto atom2 = atomic_structure[neighbour2.index];

                    auto species_triplet = SpeciesTriplet{atom0.species(),{atom1.species(), atom2.species()}};

                    if (!indexes.contains(species_triplet)) {
                        indexes[species_triplet] = ThreeBodyIndex();
                    }

                    std::array r_ij = {
                        atom1.position() + neighbour1.offset - atom0.position(),
                        atom2.position() + neighbour2.offset - atom0.position(),
                        atom2.position() + neighbour2.offset - (atom1.position() + neighbour1.offset)
                    };
                    indexes[species_triplet].push_back({
                        .atom_index = {atomIndex, neighbour1.index, neighbour2.index},
                        .r_ij = r_ij,
                        .grad_rij_wrt_rj = {r_ij[0].normalize(), r_ij[1].normalize(), r_ij[2].normalize()},
                        .fCut01 = cutoff_function_->evaluate(r_ij[0].len()),
                        .fCut02 = cutoff_function_->evaluate(r_ij[1].len()),
                        .dfCut_dr_01 = cutoff_function_->differentiate(r_ij[0].len()),
                        .dfCut_dr_02 = cutoff_function_->differentiate(r_ij[1].len()),
                        .q = toInvariantTriplet(r_ij[0].len(), r_ij[1].len(), r_ij[2].len()),
                        .dq_k_dr_ij = invariantTripletGradients(r_ij[0].len(), r_ij[1].len())
                    });
                }
            }
        }

        return indexes;
    }

    bool ThreeBodyDescriptorFinder::checkSpecies(const SpeciesTriplet& tripletInData, DataNode filters) {

        if (!filters.is_array()) {
            JGAP_LOG_AND_THROW("3b species filter is non-array: {}", filters.dump());
        }
        if (filters.empty()) return true;

        if (filters[0].is_string()) {
            // "species": ["Fe", "Ni", "Cr"] => "species": [["Fe", "Ni", "Cr"]]
            filters = DataNode::array({filters});
        }

        bool passedAFilter = false;
        for (const auto &filter: filters) {
            if (!filter.is_array() || filter.size() != 3) {
                JGAP_LOG_AND_THROW("Wrong 3b species specs: {}", filters.dump());
            }
            if (tripletInData == SpeciesTriplet(filter[0], {filter[1], filter[2]})) {
                passedAFilter = true;
                break;
            }
        }
        return passedAFilter;
    }
}
