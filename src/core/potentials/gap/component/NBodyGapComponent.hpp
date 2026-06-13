#ifndef JGAP_NBODYGAPCOMPONENT_HPP
#define JGAP_NBODYGAPCOMPONENT_HPP

#include <tbb/parallel_for.h>
#include "core/atomic/energy/AtomicQuantities.hpp"
#include "core/kernels/Kernel.hpp"
#include "core/kernels/SquaredExpKernel.hpp"
#include "core/transform/ClusterTransformation.hpp"
#include "GapComponent.hpp"
#include "core/Matrix.hpp"

namespace jgap {

    template<size_t Dim, size_t ClusterSize, ClusterSymmetry ClusterSym,
            CKernelOfDim<Dim> TKernel>
    class NBodyGapComponent : public GapComponent {
    public:
        static constexpr size_t Dependencies = Cluster<ClusterSize>::NSeparations;

        NBodyGapComponent(const SpeciesSet<ClusterSize, ClusterSym> species,
                          const ValuePtr<ClusterTransformation<Dim, ClusterSize> >& transformation,
                          const TKernel& kernel,
                          const std::vector<Descriptor<Dim> >& sparse_points,
                          const std::vector<Real>& optional_coeffs = {})
            : species(species),
              transformation(transformation),
              symmetry_factor(transformation->symmetryFactor()),
              kernel(kernel),
              sparse_points(sparse_points) {

            if (!optional_coeffs.empty()) {
                setCoefficients(optional_coeffs);
            }
        }

        std::optional<AtomicQuantities> covariate(const NeighbourList& neighbour_list) const override {

            auto clusters = neighbour_list.findAllClusters<WithDerivatives>(species);
            if (clusters.empty()) {
                return std::nullopt;
            }

            AtomicQuantities result(nSparsePoints(), neighbour_list.nAtoms());

            for (const auto& cluster: clusters) {
                auto descriptor = transformation->evaluateAndDifferentiate(cluster);

                for (size_t sparse_idx = 0; sparse_idx < sparse_points.size(); sparse_idx++) {
                    auto [K, gradK_wrt_q] = kernel.valueAndGradient(
                        sparse_points[sparse_idx].value, descriptor.value
                    );

                    K *= symmetry_factor;
                    for(auto& grad_val: gradK_wrt_q) {
                        grad_val *= symmetry_factor;
                    }

                    result.energy(sparse_idx) += K;

                    for (size_t i = 0; i < ClusterSize; i++) {
                        for (size_t j = i + 1; j < ClusterSize; j++) {
                            const auto sep_idx = flattenedIndex(i, j);
                            const auto& separation_deriv = cluster.derivativesBetween(i, j);
                            const auto& derivatives_wrt_r_norms = descriptor.derivatives[sep_idx];

                            Real dK_drij_norm = 0.0;
                            for (size_t dim = 0; dim < Dim; dim++) {
                                dK_drij_norm += derivatives_wrt_r_norms[dim] * gradK_wrt_q[dim];
                            }

                            result.force(sparse_idx, cluster.atom_indexes[i]) += separation_deriv.direction * dK_drij_norm;
                            result.force(sparse_idx, cluster.atom_indexes[j]) -= separation_deriv.direction * dK_drij_norm;

                            result.virials(sparse_idx) += separation_deriv.virials * dK_drij_norm;
                        }
                    }
                }
            }

            return result;
        }

        Matrix sparseToSparseCovariance() const override {
            Matrix result(nSparsePoints(), nSparsePoints());
            for (size_t i = 0; i < nSparsePoints(); i++) {
                for (size_t j = i; j < nSparsePoints(); j++) {
                    result(i, j) = kernel.value(sparse_points[i].value, sparse_points[j].value);
                    result(j, i) = result(i, j);
                }
            }
            return result;
        }

        size_t nSparsePoints() const override {
            return sparse_points.size();
        }

        Cutoffs getCutoffs() const override {
            return transformation->getCutoffs();
        }

        std::unique_ptr<GapComponent> clone() const override {
            return std::make_unique<NBodyGapComponent>(*this);
        }

        void tabulate(TabulationData &tables) const override {

            Grid<Dependencies>* table_ref = nullptr;

            if constexpr (ClusterSize == 2 && ClusterSym == Symmetric) {
                table_ref = &tables.two_body_grids.getValueGrid(species);
            } else if constexpr (ClusterSize == 3 && ClusterSym == HasCentralAtom) {
                table_ref = &tables.three_body_grids.getValueGrid(species);
            } else if (table_ref == nullptr) {
                JGAP_LOG_AND_THROW("Tabulation not implemented");
            }

            auto& table = *table_ref;

            tbb::parallel_for(static_cast<size_t>(0), table.data_flat.size(), [&](size_t i) {
                auto indices = table.getIndices(i);
                auto pos = table.getCoord(indices);

                auto cluster = TabulationData::gridPosAsCluster<ClusterSize>(pos);
                auto transformed = transformation->evaluate(cluster);

                Real& value = table.data_flat[i];
                for (const auto& sparse_point : sparse_points) {
                    value += kernel.value(sparse_point.value, transformed.value) * symmetry_factor;
                }
            });
        }

    private:
        SpeciesSet<ClusterSize, ClusterSym> species;
        ValuePtr<ClusterTransformation<Dim, ClusterSize>> transformation;
        Real symmetry_factor;
        TKernel kernel;
        std::vector<Descriptor<Dim>> sparse_points;
    };
}

#endif