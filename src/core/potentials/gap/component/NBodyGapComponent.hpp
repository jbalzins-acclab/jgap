#ifndef JGAP_NBODYGAPCOMPONENT_HPP
#define JGAP_NBODYGAPCOMPONENT_HPP

#include <tbb/parallel_for.h>
#include <set>
#include "core/atomic/energy/AtomicQuantities.hpp"
#include "core/kernels/Kernel.hpp"
#include "core/kernels/SquaredExpKernel.hpp"
#include "core/transform/NBodyTransformation.hpp"
#include "core/sparsification/Sparsifier.hpp"
#include "GapComponent.hpp"
#include "core/Matrix.hpp"

namespace jgap {

    template<size_t Dim, size_t ClusterSize, ClusterSymmetry ClusterSym, CKernelOfDim<Dim> TKernel>
    requires CClusterFindiningImplemented<ClusterSize, ClusterSym>
    class NBodyGapComponent : public GapComponent {
    public:
        static constexpr size_t Dependencies = Cluster<ClusterSize>::NSeparations;

        static std::vector<NBodyGapComponent> createComponents(
            const std::vector<Atoms>& training_data,
            const ValuePtr<NBodyTransformation<Dim, ClusterSize>>& transformation,
            const TKernel& kernel,
            const Sparsifier<Dim>& sparsifier) {

            std::set<SpeciesSet<ClusterSize, ClusterSym>> all_species_sets;
            Real cutoff = transformation->getCutoffs().maxOverall();

            for (const auto& atoms: training_data) {
                NeighbourList nl(atoms, cutoff);
                auto sets = nl.getSpeciesSets<ClusterSize, ClusterSym>();
                all_species_sets.insert(sets.begin(), sets.end());
            }

            std::vector<NBodyGapComponent> components;
            for (const auto& species_set : all_species_sets) {
                components.emplace_back(species_set, transformation, kernel, sparsifier, training_data);
            }

            return components;
        }

        NBodyGapComponent(const SpeciesSet<ClusterSize, ClusterSym> species,
                          const ValuePtr<NBodyTransformation<Dim, ClusterSize>>& transformation,
                          const TKernel& kernel,
                          const std::vector<Descriptor<Dim>>& sparse_points,
                          const std::vector<Real>& optional_coeffs = {})
            : species(species),
              transformation(transformation),
              symmetry_factor(transformation->symmetryFactor()),
              kernel(kernel),
              sparse_points(sparse_points) {

            if (!optional_coeffs.empty()) {
                this->setCoefficients(optional_coeffs);
            }
        }

        NBodyGapComponent(const SpeciesSet<ClusterSize, ClusterSym> species,
                          const ValuePtr<NBodyTransformation<Dim, ClusterSize> >& transformation,
                          const TKernel& kernel,
                          const Sparsifier<Dim>& sparsifier,
                          const std::vector<Atoms>& training_data,
                          const std::vector<Real>& optional_coeffs = {})
            : NBodyGapComponent(species, transformation, kernel,
                            sparsifier.selectSparsePoints(getAllDescriptors(training_data, species, transformation)),
                            optional_coeffs) {
        }

        std::optional<AtomicQuantities> covariate(const NeighbourList& neighbour_list) const override {

            auto clusters = neighbour_list.findAllClusters<WithGradients>(species);
            if (clusters.empty()) {
                return std::nullopt;
            }
            ClusterFinder::find<WithGradients>();
            AtomicQuantities result(nSparsePoints(), neighbour_list.nAtoms());

            for (const auto& cluster: clusters) {
                auto descriptor = transformation->evaluateAndDifferentiate(cluster);

                for (size_t sparse_idx = 0; sparse_idx < sparse_points.size(); sparse_idx++) {
                    auto [K, gradK_wrt_q] = kernel.valueAndGradient(
                        sparse_points[sparse_idx], descriptor
                    );


                    MOVE to THE END
                    K *= symmetry_factor;
                    for (auto& grad_val: gradK_wrt_q) {
                        grad_val *= symmetry_factor;
                    }

                    result.energy(sparse_idx) += K;

                    for (size_t i = 0; i < ClusterSize; i++) {
                        for (size_t j = i + 1; j < ClusterSize; j++) {
                            const auto sep_idx = flattenedIndex(i, j);
                            const auto& separation_derivatives = cluster.derivativesBetween(i, j);
                            const auto& derivatives_wrt_r_norms = descriptor.derivatives[sep_idx];

                            Real dK_drnorm = 0.0;
                            for (size_t dim = 0; dim < Dim; dim++) {
                                dK_drnorm += derivatives_wrt_r_norms[dim] * gradK_wrt_q[dim];
                            }

                            result.force(sparse_idx, cluster.atom_indexes[i])
                                += separation_derivatives.direction * dK_drnorm;
                            result.force(sparse_idx, cluster.atom_indexes[j])
                                -= separation_derivatives.direction * dK_drnorm;

                            result.virials(sparse_idx) += separation_derivatives.virials * dK_drnorm;
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
                    result(i, j) = kernel.value(sparse_points[i], sparse_points[j]);
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

        NBodyGapComponent* clone() const override {
            return new NBodyGapComponent(*this);
        }

        const SpeciesSet<ClusterSize, ClusterSym>& getSpecies() const { return species; }
        const ValuePtr<NBodyTransformation<Dim, ClusterSize>>& getTransformation() const { return transformation; }
        const TKernel& getKernel() const { return kernel; }
        const std::vector<Descriptor<Dim>>& getSparsePoints() const { return sparse_points; }

        void tabulate(TabulationData &tables) const override {
            const auto& coeffs = this->getCoefficients();
            if (coeffs.empty()) {
                JGAP_LOG_AND_THROW("Coefficients must be set before tabulation");
            }

            Grid<Dependencies>* table_ref = nullptr;

            if constexpr (ClusterSize == 2 && ClusterSym == FullSymmetry) {
                table_ref = &tables.two_body_grids.getValueGrid(species);
            } else if constexpr (ClusterSize == 3 && ClusterSym == NodeSymmetric) {
                table_ref = &tables.three_body_grids.getValueGrid(species);
            } else {
                JGAP_LOG_AND_THROW("Tabulation not implemented");
            }

            auto& table = *table_ref;

            tbb::parallel_for(static_cast<size_t>(0), table.data_flat.size(), [&](size_t i) {
                auto indices = table.getIndices(i);
                auto pos = table.getCoord(indices);

                auto cluster = TabulationData::gridPosAsCluster<ClusterSize>(pos);
                auto transformed = transformation->evaluate(cluster);

                Real& value = table.data_flat[i];
                for (size_t sparse_idx = 0; sparse_idx < sparse_points.size(); sparse_idx++) {
                    value += coeffs[sparse_idx]
                                * kernel.value(sparse_points[sparse_idx], transformed)
                                * symmetry_factor;
                }
            });
        }

    private:
        static std::vector<Descriptor<Dim>> getAllDescriptors(
            const std::vector<Atoms>& training_data,
            const SpeciesSet<ClusterSize, ClusterSym>& species,
            const ValuePtr<NBodyTransformation<Dim, ClusterSize>>& transformation)
        {
            std::vector<Descriptor<Dim>> all_descriptors;
            for (const auto& atoms : training_data) {
                NeighbourList nl(atoms, transformation->getCutoffs().maxOverall());
                auto clusters = nl.findAllClusters<ValueOnly>(species);
                for (const auto& cluster: clusters) {
                    all_descriptors.push_back(transformation->evaluate(cluster));
                }
            }
            return all_descriptors;
        }

        SpeciesSet<ClusterSize, ClusterSym> species;
        ValuePtr<NBodyTransformation<Dim, ClusterSize>> transformation;
        Real symmetry_factor;
        TKernel kernel;
        std::vector<Descriptor<Dim>> sparse_points;
    };
}

#endif