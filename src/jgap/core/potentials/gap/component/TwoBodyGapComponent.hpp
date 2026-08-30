#ifndef JGAP_TWOBODYGAPCOMPONENT_HPP
#define JGAP_TWOBODYGAPCOMPONENT_HPP

#include <cassert>
#include <ranges>
#include <set>
#include "GapComponent.hpp"
#include "jgap/core/Matrix.hpp"
#include "jgap/core/ValuePtr.hpp"
#include "jgap/core/atomic/energy/AtomicQuantities.hpp"
#include "jgap/core/atomic/iteration/Cluster2Expansion.hpp"
#include "jgap/core/atomic/neighbours/NeighbourLists.hpp"
#include "jgap/core/kernels/Kernel.hpp"
#include "jgap/core/sparsification/Sparsifier.hpp"
#include "jgap/core/transform/nbody/2b/TwoBodyTransformation.hpp"
#include "jgap/utils/Utils.hpp"

namespace jgap {

    template<size_t Dim, CKernelOfDim<Dim> TKernel>
    class TwoBodyGapComponent : public GapComponent {
    public:
        static constexpr size_t Dependencies = 1;

        TwoBodyGapComponent(
            const Species2Sorted& species,
            ValuePtr<TwoBodyTransformation<Dim>> transformation,
            TKernel kernel,
            std::vector<Descriptor<Dim>> sparse_points,
            const std::vector<Real>& optional_coeffs = {}
        ) :
            species(species),
            transformation(std::move(transformation)),
            kernel(std::move(kernel)),
            sparse_points(std::move(sparse_points)) {
            if (!optional_coeffs.empty()) {
                this->setCoefficients(optional_coeffs);
            }
        }

        TwoBodyGapComponent(
            const Species2Sorted species,
            const ValuePtr<TwoBodyTransformation<Dim>>& transformation,
            const TKernel& kernel,
            const Sparsifier<Dim>& sparsifier,
            const std::vector<Atoms>& training_data
        ) :
            species(species), transformation(transformation), kernel(kernel) {
            auto all_descriptors = getAllDescriptors(training_data, species, transformation);
            sparse_points = sparsifier.selectSparsePoints(all_descriptors);
        }

        std::optional<AtomicQuantities> covariate(const NeighbourLists& neighbour_list) const override {
            AtomicQuantities result(nSparsePoints(), neighbour_list.nAtoms());

            Cluster2Expansion expansion(species);

            bool found = expansion.forEach(neighbour_list, [&](const Cluster2& cluster) {
                auto descriptor = transformation->evaluateAndDifferentiate(cluster);

                for (size_t sparse_idx = 0; sparse_idx < sparse_points.size(); sparse_idx++) {
                    auto [K, gradK_wrt_q] = kernel.valueAndGradient(sparse_points[sparse_idx], descriptor.value);

                    result.energy(sparse_idx) += K;

                    Vector3 f1{};
                    for (size_t dim = 0; dim < Dim; dim++) {
                        f1 -= gradK_wrt_q[dim] * descriptor.grad_r1[dim];
                    }

                    result.force(sparse_idx, cluster.idx1) += f1;
                    result.force(sparse_idx, cluster.idx0) -= f1;

                    result.virials(sparse_idx) += Virials::dyadic(cluster.separation01.vec(), f1);
                }
            });

            if (!found) {
                return std::nullopt;
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
            return {species.nodes[0], species.nodes[1]};
        }

        TwoBodyGapComponent* clone() const override { return new TwoBodyGapComponent(*this); }

        const Species2Sorted& getSpecies() const { return species; }
        const ValuePtr<TwoBodyTransformation<Dim>>& getTransformation() const { return transformation; }
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

            Grid<Dependencies>* table_ref = &tables.two_body_grids.getValueGrid(species);
            auto& table = *table_ref;

            for (size_t i = 0; i < table.data_flat.size(); ++i) {
                auto indices = table.getIndices(i);
                auto pos = table.getCoord(indices);

                auto cluster = TabulationData::gridPosAsCluster2(pos);
                auto transformed = transformation->evaluate(cluster);

                Real& value = table.data_flat[i];
                for (size_t sparse_idx = 0; sparse_idx < sparse_points.size(); sparse_idx++) {
                    constexpr Real iteration_reduction_factor = 2.0;
                    value += coeffs[sparse_idx] * kernel.value(sparse_points[sparse_idx], transformed)
                             * iteration_reduction_factor;
                }
            }
        }

    private:
        static std::vector<Descriptor<Dim>> getAllDescriptors(
            const std::vector<Atoms>& training_data,
            const Species2Sorted& species,
            const ValuePtr<TwoBodyTransformation<Dim>>& transformation
        ) {
            std::vector<Descriptor<Dim>> all_descriptors;
            Cluster2Expansion exp(species);

            for (const auto& atoms: training_data) {
                NeighbourLists nl(atoms, transformation->getCutoffs().maxOverall());
                exp.forEach(nl, [&](const Cluster2& cluster) {
                    all_descriptors.push_back(transformation->evaluate(cluster));
                });
            }
            return all_descriptors;
        }

        Species2Sorted species;
        ValuePtr<TwoBodyTransformation<Dim>> transformation;
        TKernel kernel;
        std::vector<Descriptor<Dim>> sparse_points;
    };
}

#endif
