#include "core/descriptors/3b/ThreeBodyDescriptorFinder.hpp"

#include <random>
#include <tbb/parallel_for_each.h>

#include "io/log/StdoutLogger.hpp"
#include "utils/Utils.hpp"

namespace jgap {

    /*
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
                    sqrt(std::max/*numeric safety*//*(
                        pow(iter.pos.x, 2) + pow(iter.pos.y, 2) - 2.0 * iter.pos.x * iter.pos.y * iter.pos.z, 0.0
                        ))
                );

                double contribution = 0.0;
                for (const size_t kernelId: kernelIds) {
                    contribution += _kernels[kernelId]->value(
                        {.q = invariantTriplet, .f_cut = invariantTripletToCutoff(invariantTriplet)}
                        ) * 2.0/*q_ijk + q_jik*//* * _kernels[kernelId]->coefficient.value();
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
    */

    ThreeBodyDescriptorsFiltered ThreeBodyDescriptorFinder::findSeparations(
        const AtomicStructure &atomic_structure) const {

        const double cutoff = cutoff_function->getCutoff();

        std::map<EncodedSpeciesSets, std::vector<NewThreeBodyDescriptor>> result;

        for (size_t atom0_index = 0; atom0_index < atomic_structure.size(); atom0_index++) {
            auto atom0 = atomic_structure[atom0_index];

            for (size_t nl_index1 = 0; nl_index1 < atom0.neighboursAscendingSeparation().size(); nl_index1++) {

                auto neighbour1 = atom0.neighboursAscendingSeparation()[nl_index1];
                auto atom1 = atomic_structure[neighbour1.index];
                if (neighbour1.distance > cutoff) continue;

                for (size_t nl_index2 = nl_index1 + 1; nl_index2 < atom0.neighboursAscendingSeparation().size(); nl_index2++) {

                    auto neighbour2 = atom0.neighboursAscendingSeparation()[nl_index2];
                    auto atom2 = atomic_structure[neighbour2.index];
                    if (neighbour2.distance > cutoff) continue;

                    auto species_triplet = SpeciesEncoder::asSet(
                        atom0.speciesEncoded(),
                        atom1.speciesEncoded(), atom2.speciesEncoded()
                        );

                    if (!result.contains(species_triplet)) {
                        result[species_triplet] = {};
                    }

                    std::array r_ij = {
                        atom1.position() + neighbour1.offset - atom0.position(), // 0->1
                        atom2.position() + neighbour2.offset - atom0.position(), // 0->2
                        atom2.position() + neighbour2.offset - (atom1.position() + neighbour1.offset) // 1->2
                    };

                    double r_01 = r_ij[0].len();
                    double r_02 = r_ij[1].len();
                    double r_12 = r_ij[1].len();

                    double f_cut_01 = cutoff_function->evaluate(r_01);
                    double f_cut_02 = cutoff_function->evaluate(r_02);

                    if (!calculate_derivatives) {
                        result[species_triplet].push_back(NewThreeBodyDescriptor{
                            .value = { r_01, r_02, r_12 },
                            .f_cut = f_cut_01 * f_cut_02
                        });
                        continue;
                    }
                    // ELSE

                    std::array grad_rij_wrt_j/*[j][ij]*/ = {
                        // ij = 01, 02, 12
                        std::array{ r_ij[0].normalize() * -1.0, r_ij[1].normalize() * -1.0, Vector3() },
                        std::array{ r_ij[0].normalize(),  Vector3(), r_ij[2].normalize() * -1.0 },
                        std::array{ Vector3(), r_ij[1].normalize(), r_ij[2].normalize() }
                    };

                    std::array<Virials, 3> virials{};

                    #pragma unroll
                    for (size_t ij = 0; ij < 3; ij++) {
                        constexpr std::array i_ij = {0, 0, 1};
                        virials[ij] = {
                            .xx = r_ij[ij].x * grad_rij_wrt_j[i_ij[ij]][ij].x,
                            .xy = r_ij[ij].x * grad_rij_wrt_j[i_ij[ij]][ij].y,
                            .xz = r_ij[ij].x * grad_rij_wrt_j[i_ij[ij]][ij].z,
                            .yy = r_ij[ij].y * grad_rij_wrt_j[i_ij[ij]][ij].y,
                            .yz = r_ij[ij].y * grad_rij_wrt_j[i_ij[ij]][ij].z,
                            .zz = r_ij[ij].z * grad_rij_wrt_j[i_ij[ij]][ij].z,
                        };
                    }

                    double df_cut_01_dr = cutoff_function->differentiate(r_01);
                    double df_cut_02_dr = cutoff_function->differentiate(r_02);

                    Vector3 df_cut_wrt1  = grad_rij_wrt_j[1][0] * df_cut_01_dr * f_cut_02;
                    Vector3 df_cut_wrt2  = grad_rij_wrt_j[2][1] * df_cut_02_dr * f_cut_01;
                    Vector3 df_cut_wrt0 = (df_cut_wrt1 + df_cut_wrt2) * -1;

                    result[species_triplet].push_back(NewThreeBodyDescriptor{
                        .value = { r_01, r_02, r_12 },
                        .virials = virials,
                        .gradients = {
                            NewThreeBodyDescriptor::GradientData{ atom0_index, grad_rij_wrt_j[0] },
                            NewThreeBodyDescriptor::GradientData{ neighbour1.index, grad_rij_wrt_j[1] },
                            NewThreeBodyDescriptor::GradientData{ neighbour2.index, grad_rij_wrt_j[2] }
                        },
                        .f_cut = f_cut_01 * f_cut_02,
                        .f_cut_gradients = { df_cut_wrt0, df_cut_wrt1, df_cut_wrt2 }
                    });
                }
            }
        }

        return result;
    }

    ThreeBodyDescriptorsFiltered ThreeBodyDescriptorFinder::find(const AtomicStructure &atomic_structure) const {
        const auto separations = findSeparations(atomic_structure);
        return toInvariantTriplets(separations);
    }

    ThreeBodyDescriptorsFiltered ThreeBodyDescriptorFinder::toInvariantTriplets(
        const ThreeBodyDescriptorsFiltered& separations) const {

        ThreeBodyDescriptorsFiltered result{};
        for (const auto& [filter, per_filter]: separations) {
            result[filter] = {};
            for (const auto& descriptor: per_filter) {
                result[filter].push_back(toInvariantTriplet(descriptor));
            }
        }
        return result;
    }

    NewThreeBodyDescriptor ThreeBodyDescriptorFinder::toInvariantTriplet(
        const NewThreeBodyDescriptor &separations) const {

        const auto& [q, dq_dr_ij] = toInvariantTriplet(separations.value);

        if (!calculate_derivatives) {
            return {
                .value = q,
                .f_cut = separations.f_cut
            };
        }
        // ELSE

        std::array<Virials, 3> virials{}; // [k]
        for (size_t ij = 0; ij < 3; ij++) {
            virials[0] += separations.virials[0] * dq_dr_ij[ij].x;
            virials[1] += separations.virials[1] * dq_dr_ij[ij].y;
            virials[2] += separations.virials[2] * dq_dr_ij[ij].z;
        }

        std::array<NewThreeBodyDescriptor::GradientData, 3>  dqk_drl{}; // [l].grad(wrt l)[k]
        for (size_t l = 0; l < 3; l++) {
            auto& dq = dqk_drl[l];
            dq.wrt_atom_index = separations.gradients[l].wrt_atom_index;

            const auto& drij_drl = separations.gradients[l].gradients;

            auto& g0 = dq.gradients[0];
            auto& g1 = dq.gradients[1];
            auto& g2 = dq.gradients[2];

            for (size_t ij = 0; ij < 3; ij++) {
                g0 += drij_drl[ij] * dq_dr_ij[ij].x;
                g1 += drij_drl[ij] * dq_dr_ij[ij].y;
                g2 += drij_drl[ij] * dq_dr_ij[ij].z;
            }
        }

        return {
            .value = q,
            .virials = virials,
            .gradients = dqk_drl,
            .f_cut = separations.f_cut,
            .f_cut_gradients = separations.f_cut_gradients
        };
    }

    ThreeBodySeparationsTransformed ThreeBodyDescriptorFinder::toInvariantTriplet(
        const std::array<double, 3>& separations) const {

        const auto&[r01, r02, r12] = separations;

        std::array q = { r01 + r02, (r01-r02) * (r01 - r02), r12 }; // [k]

        if (!calculate_derivatives) return {.q = q};

        std::array dq_dr_ij = { // [ij]
            Vector3{1.0, 2.0 * (r01 - r02), 0.0},
            Vector3{1.0, 2.0 * (r02 - r01), 0.0},
            Vector3{0.0, 0.0,               1.0}
        };

        return { q, dq_dr_ij };
    }
}
