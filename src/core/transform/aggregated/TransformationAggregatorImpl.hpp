#ifndef JGAP_TRANSFORMATIONAGGREGATORIMPL_HPP
#define JGAP_TRANSFORMATIONAGGREGATORIMPL_HPP
#include "TransformationAggregator.hpp"
#include "core/tabulation/TabulationData.hpp"

namespace jgap {

    template<size_t Dim, size_t ClusterSize>
    class TransformationAggregatorImpl final : public TransformationAggregator<Dim> {
    public:
        using TransformationPtr = ValuePtr<ClusterTransformation<Dim, ClusterSize>>;

        TransformationAggregatorImpl(const Species central_atom_species)
            : central_atom_species(central_atom_species) {
        }

        void extend(SpeciesSet<ClusterSize, HasCentralAtom> species_set, TransformationPtr transformation) {

            if (species_set.getRoot() != central_atom_species) {
                JGAP_LOG_AND_THROW("Transformation intended for a cluster "
                                   "whose root doesnt match central atom species");
            }

            transformations.insert({species_set, std::move(transformation)});
        }

        std::map<size_t, ManyBodyDescriptor<Dim>> aggregate(const NeighbourList& nl) const override {
            std::map<size_t, ManyBodyDescriptor<Dim>> aggregated_descriptors;

            if (!nl.atoms_by_species.contains(central_atom_species)) {
                return aggregated_descriptors;
            }

            const auto& central_atom_indices = nl.atoms_by_species.at(central_atom_species);

            for (size_t idx: central_atom_indices) {
                aggregated_descriptors[idx] = ManyBodyDescriptor<Dim>(nl.nAtoms());
            }

            // Stage 1: Form full descriptors per central atom
            for (const auto& [species_set, transformation] : transformations) {
                auto clusters = nl.findAllClusters<WithDerivatives>(species_set);
                for (const auto& cluster: clusters) {
                    size_t central_idx = cluster.atom_indexes[0];

                    auto& descriptor = aggregated_descriptors.at(central_idx);
                    auto contribution = transformation->evaluateAndDifferentiate(cluster);

                    // Accumulate value
                    for (size_t d = 0; d < Dim; d++) {
                        descriptor.value[d] += contribution.value[d];
                    }

                    // Accumulate forces and virials directly via chain rule wrt internal variables
                    for (size_t i = 0; i < ClusterSize; i++) {
                        for (size_t j = i + 1; j < ClusterSize; j++) {
                            const auto sep_idx = flattenedIndex(i, j);
                            const auto& separation = cluster.derivativesBetween(i, j);
                            const auto& derivs = contribution.derivatives[sep_idx];

                            for (size_t dim = 0; dim < Dim; dim++) {
                                Vector3 force_contrib = separation.direction * derivs[dim];
                                // Force is the negative gradient: F_i = -dE/dr_i = dE/dr_ij * direction
                                descriptor.forces[cluster.atom_indexes[i]][dim] += force_contrib;
                                // F_j = -dE/dr_j = -dE/dr_ij * direction
                                descriptor.forces[cluster.atom_indexes[j]][dim] -= force_contrib;

                                descriptor.virials[dim] += separation.virials * derivs[dim];
                            }
                        }
                    }
                }
            }
            return aggregated_descriptors;
        }

        Cutoffs getCutoffs() const override {
            Cutoffs combined;
            for (const auto& trans : transformations | std::views::values) {
                combined += trans->getCutoffs();
            }
            return combined;
        }

        std::unique_ptr<TransformationAggregator<Dim>> clone() const override {
            return std::make_unique<TransformationAggregatorImpl>(*this);
        }

        void tabulateNewManyBodyGrid(TabulationData &tables) const override {

            if constexpr (Dim == 1 && ClusterSize == 2) {

                ManyBodyGrids<2>& eam_grids = tables.newEamGrid(central_atom_species);

                for (const auto& [species_set, transformation] : transformations) {
                    auto& grid = eam_grids.aggregator_grids.getValueGrid(species_set);

                    for (auto cell: grid) {
                        Cluster<2> as_cluster{};
                        as_cluster.between(0, 1) = cell.pos[0];
                        cell.value += transformation->evaluate(as_cluster).value[0];
                    }
                }

            } else {
                JGAP_LOG_AND_THROW("Tabulation not implemented");
            }
        }

        Species getCentralSpecies() const override {
            return central_atom_species;
        }

        std::set<Species> getAllSpecies() const override {
            std::set result{central_atom_species};

            for (const auto& species_set: transformations | std::views::keys) {
                for (Species species: species_set.getNodes()) {
                    result.insert(species);
                }
            }

            return result;
        }

        const auto& getTransformations() const {
            return transformations;
        }

    private:
        Species central_atom_species;
        std::multimap<SpeciesSet<ClusterSize, HasCentralAtom>, TransformationPtr> transformations;
    };
}

#endif
