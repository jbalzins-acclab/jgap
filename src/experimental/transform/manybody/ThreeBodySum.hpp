#ifndef JGAP_THREEBODYSUM_HPP
#define JGAP_THREEBODYSUM_HPP

#include <cassert>
#include <map>
#include <ranges>
#include <set>

#include "core/atomic/descriptor/ManyBodyDescriptors.hpp"
#include "core/atomic/iteration/AtomicCluster3Expansion.hpp"
#include "core/atomic/species/composition/Species3AtomicSorted.hpp"
#include "core/tabulation/TabulationData.hpp"
#include "core/transform/manybody/NBodyAggregator.hpp"
#include "core/transform/nbody/3b/ThreeBodyTransformation.hpp"

namespace jgap {

    template<size_t Dim>
    class ThreeBodySum final : public NBodyAggregator<Dim> {
    public:
        using TransformationPtr = ValuePtr<ThreeBodyTransformation<Dim>>;

        ThreeBodySum(const Species central_atom_species) : central_atom_species(central_atom_species) {}

        void extend(Species3AtomicSorted species_set, TransformationPtr transformation) {
            if (species_set.root != central_atom_species) {
                JGAP_LOG_AND_THROW(
                    "Transformation intended for a cluster "
                    "whose root doesn't match central atom species");
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
                    AtomicCluster3Expansion expansion(species_set);
                    auto expansion_result = expansion.find(atom_index, nl, CalculationType::WithGradients);
                    assert(expansion_result.derivatives.has_value());

                    const bool needs_permutation = !transformation->isSwapInvariant(1, 2);

                    for (const auto& [cluster, cluster_derivs]:
                         std::views::zip(expansion_result.clusters, *expansion_result.derivatives)) {
                        auto contribution = transformation->evaluateAndDifferentiate(cluster);
                        aggregated_descriptors.add(descriptor_index, cluster, cluster_derivs, contribution);

                        if (needs_permutation) {
                            Cluster3 permuted_cluster = cluster;
                            std::swap(permuted_cluster.atom_indexes[1], permuted_cluster.atom_indexes[2]);
                            std::swap(permuted_cluster.separation_magnitudes[0],
                                      permuted_cluster.separation_magnitudes[1]);

                            Cluster3Derivatives permuted_derivs = cluster_derivs;
                            std::swap(permuted_derivs.val[0], permuted_derivs.val[1]);
                            permuted_derivs.val[2].direction *= -1.0;

                            auto contribution2 = transformation->evaluateAndDifferentiate(permuted_cluster);
                            aggregated_descriptors.add(descriptor_index, permuted_cluster, permuted_derivs,
                                                       contribution2);
                        }
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

        ThreeBodySum* clone() const override { return new ThreeBodySum(*this); }

        void tabulateNewManyBodyGrid(TabulationData& tables) const override {
            JGAP_LOG_AND_THROW("Tabulation not implemented for ThreeBodySum");
        }

        Species getCentralSpecies() const override { return central_atom_species; }

        std::set<Species> getAllSpecies() const override {
            std::set result{central_atom_species};
            for (const auto& species_set: transformations | std::views::keys) {
                result.insert(species_set.nodes[0]);
                result.insert(species_set.nodes[1]);
            }
            return result;
        }

        const auto& getTransformations() const { return transformations; }

    private:
        Species central_atom_species;
        std::multimap<Species3AtomicSorted, TransformationPtr> transformations;
    };
}

#endif
