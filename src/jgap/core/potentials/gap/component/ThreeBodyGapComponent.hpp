#ifndef JGAP_THREEBODYGAPCOMPONENT_HPP
#define JGAP_THREEBODYGAPCOMPONENT_HPP

#include <cassert>
#include <ranges>
#include <set>
#include "../../../transform/nbody/3b/ThreeBodyTransformation.hpp"
#include "GapComponent.hpp"
#include "jgap/core/Matrix.hpp"
#include "jgap/core/atomic/energy/AtomicQuantities.hpp"
#include "jgap/core/atomic/iteration/Cluster3Expansion.hpp"
#include "jgap/core/atomic/iteration/ClusterPermutationMode.hpp"
#include "jgap/core/kernels/Kernel.hpp"
#include "jgap/core/sparsification/Sparsifier.hpp"
#include "jgap/utils/Utils.hpp"

namespace jgap {

    template<size_t Dim, CKernelOfDim<Dim> TKernel>
    class ThreeBodyGapComponent : public GapComponent {
    public:
        static constexpr size_t Dependencies = 3;

        ThreeBodyGapComponent(
            const Species3AtomicSorted species,
            const ValuePtr<ThreeBodyTransformation<Dim>>& transformation,
            const TKernel& kernel,
            const std::vector<Descriptor<Dim>>& sparse_points,
            const std::vector<Real>& optional_coeffs = {}
        ) :
            species(species),
            transformation(transformation),
            kernel(kernel),
            sparse_points(sparse_points),
            expansion(
                species,
                transformation->isSwapInvariant(1, 2) ? ClusterPermutationMode::NoNodePermutation
                                                      : ClusterPermutationMode::PermuteSameSpeciesNodes
            ) {
            if (!optional_coeffs.empty()) {
                this->setCoefficients(optional_coeffs);
            }
        }

        ThreeBodyGapComponent(
            const Species3AtomicSorted species,
            const ValuePtr<ThreeBodyTransformation<Dim>>& transformation,
            const TKernel& kernel,
            const Sparsifier<Dim>& sparsifier,
            const std::vector<Atoms>& training_data
        ) :
            ThreeBodyGapComponent(
                species,
                transformation,
                kernel,
                sparsifier.selectSparsePoints(getAllDescriptors(training_data, species, transformation))
            ) {}

        std::optional<AtomicQuantities> covariate(const NeighbourLists& neighbour_list) const override {
            AtomicQuantities result(nSparsePoints(), neighbour_list.nAtoms());
            bool found_any = false;

            auto it = neighbour_list.atoms_by_species.find(species.root);
            if (it != neighbour_list.atoms_by_species.end()) {
                for (size_t atom_index: it->second) {
                    found_any |= covariateAtom(atom_index, neighbour_list, result);
                }
            }

            if (!found_any) {
                return std::nullopt;
            }

            const Real factor = expansion.getPermutationReductionFactor();
            if (factor != 1.0) {
                result *= factor;
            }

            return result;
        }

        Matrix<RowMajor> sparseToSparseCovariance() const override {
            Matrix<RowMajor> result(nSparsePoints(), nSparsePoints());
            for (size_t i = 0; i < nSparsePoints(); i++) {
                for (size_t j = i; j < nSparsePoints(); j++) {
                    result(i, j) = kernel.value(sparse_points[i], sparse_points[j]);
                    result(j, i) = result(i, j);
                }
            }
            return result;
        }

        size_t nSparsePoints() const override { return sparse_points.size(); }

        Cutoffs getCutoffs() const override { return transformation->getCutoffs(); }

        std::set<Species> nonZeroCovarianceFor() const override {
            return {species.root, species.nodes[0], species.nodes[1]};
        }

        ThreeBodyGapComponent* clone() const override { return new ThreeBodyGapComponent(*this); }

        const Species3AtomicSorted& getSpecies() const { return species; }
        const ValuePtr<ThreeBodyTransformation<Dim>>& getTransformation() const { return transformation; }
        const TKernel& getKernel() const { return kernel; }
        const std::vector<Descriptor<Dim>>& getSparsePoints() const { return sparse_points; }

        void tabulate(TabulationData& tables) const override {
            const auto& coeffs = this->getCoefficients();
            if (coeffs.empty()) {
                JGAP_LOG_AND_THROW("Coefficients must be set before tabulation");
            }
            if (!transformation->isRotationallyInvariant()) {
                JGAP_LOG_AND_THROW("Transformation is not rotationally invariant and cannot be tabulated");
            }

            auto& table = tables.three_body_grids.getValueGrid(species);
            const Real iteration_reduction_factor = expansion.getPermutationReductionFactor();

            for (size_t i = 0; i < table.data_flat.size(); ++i) {
                auto indices = table.getIndices(i);
                auto pos = table.getCoord(indices);

                auto cluster = TabulationData::gridPosAsCluster3(pos);
                auto transformed = transformation->evaluate(cluster);

                Real& value = table.data_flat[i];
                for (size_t sparse_idx = 0; sparse_idx < sparse_points.size(); sparse_idx++) {
                    value += coeffs[sparse_idx] * kernel.value(sparse_points[sparse_idx], transformed)
                             * iteration_reduction_factor;
                }
            }
        }

    private:
        bool covariateAtom(size_t atom_index, const NeighbourLists& neighbour_list, AtomicQuantities& result) const {
            return expansion.forEach(atom_index, neighbour_list, [&](const Cluster3& cluster) {
                auto descriptor = transformation->evaluateAndDifferentiate(cluster);
                const auto r01 = cluster.separation01.vec();
                const auto r02 = cluster.separation02.vec();

                for (size_t sparse_idx = 0; sparse_idx < sparse_points.size(); sparse_idx++) {
                    auto [K, gradK_wrt_q] = kernel.valueAndGradient(sparse_points[sparse_idx], descriptor.value);

                    result.energy(sparse_idx) += K;

                    Vector3 f1{0.0_r, 0.0_r, 0.0_r};
                    Vector3 f2{0.0_r, 0.0_r, 0.0_r};
                    for (size_t dim = 0; dim < Dim; dim++) {
                        f1 -= gradK_wrt_q[dim] * descriptor.grad_r1[dim];
                        f2 -= gradK_wrt_q[dim] * descriptor.grad_r2[dim];
                    }

                    result.force(sparse_idx, cluster.idx1) += f1;
                    result.force(sparse_idx, cluster.idx2) += f2;
                    result.force(sparse_idx, cluster.idx0) -= (f1 + f2);

                    result.virials(sparse_idx) += Virials::dyadic(r01, f1);
                    result.virials(sparse_idx) += Virials::dyadic(r02, f2);
                }
            });
        }

        static std::vector<Descriptor<Dim>> getAllDescriptors(
            const std::vector<Atoms>& training_data,
            const Species3AtomicSorted& species,
            const ValuePtr<ThreeBodyTransformation<Dim>>& transformation
        ) {
            std::vector<Descriptor<Dim>> all_descriptors;
            const auto mode = transformation->isSwapInvariant(1, 2) ? ClusterPermutationMode::NoNodePermutation
                                                                    : ClusterPermutationMode::PermuteSameSpeciesNodes;
            Cluster3Expansion exp(species, mode);
            for (const auto& atoms: training_data) {
                NeighbourLists nl(atoms, transformation->getCutoffs().maxOverall());

                exp.forEach(nl, [&](const Cluster3& cluster) {
                    all_descriptors.push_back(transformation->evaluate(cluster));
                });
            }
            return all_descriptors;
        }

        Species3AtomicSorted species;
        ValuePtr<ThreeBodyTransformation<Dim>> transformation;
        TKernel kernel;
        std::vector<Descriptor<Dim>> sparse_points;
        Cluster3Expansion expansion;
    };
}

#endif
