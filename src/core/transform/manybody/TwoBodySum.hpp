#ifndef JGAP_TWOBODYSUM_HPP
#define JGAP_TWOBODYSUM_HPP

#include <cassert>
#include <map>
#include <ranges>
#include <set>

#include "NBodyAggregator.hpp"
#include "core/atomic/descriptor/ManyBodyDescriptors.hpp"
#include "core/atomic/iteration/AtomicCluster2Expansion.hpp"
#include "core/atomic/species/composition/Species2Atomic.hpp"
#include "core/tabulation/TabulationData.hpp"
#include "core/transform/nbody/2b/TwoBodyTransformation.hpp"

namespace jgap {

    template<size_t Dim>
    class TwoBodySum final : public NBodyAggregator<Dim> {
    public:
        using TransformationPtr = ValuePtr<TwoBodyTransformation<Dim>>;

        TwoBodySum(const Species central_atom_species) : central_atom_species(central_atom_species) {}

        void extend(Species2Atomic species_set, TransformationPtr transformation) {
            if (species_set.root != central_atom_species) {
                JGAP_LOG_AND_THROW(
                    "Transformation intended for a cluster "
                    "whose root doesnt match central atom species");
            }
            transformations.insert({species_set, std::move(transformation)});
        }

        ManyBodyDescriptors<Dim> aggregate(const NeighbourLists& nl) const override {
            if (!nl.atoms_by_species.contains(central_atom_species)) {
                return {};
            }
            auto& atom_indexes = nl.atoms_by_species.at(central_atom_species);

            ManyBodyDescriptors<Dim> aggregated_descriptors(atom_indexes.size(), nl.nAtoms());
            for (size_t descriptor_index = 0; descriptor_index < atom_indexes.size(); ++descriptor_index) {
                size_t atom_index = atom_indexes[descriptor_index];
                for (const auto& [species_set, transformation]: transformations) {
                    AtomicCluster2Expansion expansion(species_set);
                    auto expansion_result = expansion.find(atom_index, nl, CalculationType::WithGradients);
                    assert(expansion_result.derivatives.has_value());

                    for (const auto& [cluster, cluster_derivs]:
                         std::views::zip(expansion_result.clusters, *expansion_result.derivatives)) {
                        auto contribution = transformation->evaluateAndDifferentiate(cluster);

                        aggregated_descriptors.add(descriptor_index, cluster, cluster_derivs, contribution);
                    }
                }
            }

            return aggregated_descriptors;
        }

        Cutoffs getCutoffs() const override {
            Cutoffs combined;
            for (const auto& trans: transformations | std::views::values) {
                combined += trans->getCutoffs();
            }
            return combined;
        }

        TwoBodySum* clone() const override { return new TwoBodySum(*this); }

        void tabulateNewManyBodyGrid(TabulationData& tables) const override {
            if constexpr (Dim == 1) {
                ManyBodyGrids2<1, 1>& eam_grids = tables.newEamGrid(central_atom_species);

                for (const auto& [species_set, transformation]: transformations) {
                    auto& grid = eam_grids.aggregator_grids.getValueGrid(species_set);

                    for (auto cell: grid) {
                        Cluster2 as_cluster{};
                        as_cluster.r01 = cell.pos[0];
                        cell.value += transformation->evaluate(as_cluster)[0];
                    }
                }
            } else {
                JGAP_LOG_AND_THROW("Tabulation not implemented");
            }
        }

        Species getCentralSpecies() const override { return central_atom_species; }

        std::set<Species> getAllSpecies() const override {
            std::set result{central_atom_species};
            for (const auto& species_set: transformations | std::views::keys) {
                result.insert(species_set.node);
            }
            return result;
        }

        const auto& getTransformations() const { return transformations; }

    private:
        Species central_atom_species;
        std::multimap<Species2Atomic, TransformationPtr> transformations;
    };
}

#endif
