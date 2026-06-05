#ifndef JGAP_EXTENSIVETRANSFORMER_HPP
#define JGAP_EXTENSIVETRANSFORMER_HPP

#include "../ClusterTransformation.hpp"
#include <map>
#include <vector>

#include "core/transform/eam/EamTransformer.hpp"

namespace jgap {

    template<typename TTransformer>
    requires CClusterTransformation<TTransformer>
    class ExtensiveTransformerGroup {
    public:
        static constexpr size_t Dim = TTransformer::Dim;
        static constexpr size_t Dependencies = TTransformer::Dependencies;
        static constexpr size_t ClusterSize = TTransformer::ClusterSize;

        ExtensiveTransformerGroup(const std::map<SpeciesSet, TTransformer>& transformers)
            : transformers(transformers) {
        }

        template<CalculationType CalcType>
        auto transform(const NeighbourList &neighbour_list) const {

            double cutoff = this->getCutoff();

            TransformerGroupReturns<Dim, DYNAMIC_DEPENDENCIES, CalcType> result;

            size_t n_atoms = neighbour_list.nAtoms();
            std::vector<size_t> indexes_in_results(n_atoms);

            for (auto& [species, indices]: neighbour_list.atoms_by_species) {

                if constexpr (CalcType == CalculationType::ValueOnly) {
                    result[SpeciesSet(species)] = std::vector(
                        indices.size(), Descriptor<Dim, DYNAMIC_DEPENDENCIES, CalcType>{}
                        );
                } else {
                    result[SpeciesSet(species)] = std::vector(
                        indices.size(), Descriptor<Dim, DYNAMIC_DEPENDENCIES, CalcType>(n_atoms)
                        );
                }

                for (size_t i = 0; i < indices.size(); i++) {
                    indexes_in_results[indices[i]] = i;
                }
            }

            for (const auto& species_set: neighbour_list.getSpeciesSets<ClusterSize>()) {

                auto& result_per_species = result[species_set.rootOnly()];

                auto& transformer = transformerForSpecies(species_set);

                for (auto& cluster: neighbour_list.findAllClusters<ClusterSize>(species_set, cutoff)) {
                    auto& add_result_to = result_per_species[indexes_in_results[cluster.atom_indexes[0]]];
                    transformer.template transform<CalcType>(add_result_to, cluster);
                }
            }

            return result;
        }

    private:
        std::map<SpeciesSet, TTransformer> transformers;

        const TTransformer& transformerForSpecies(SpeciesSet species_set) const {
            return transformers.at(species_set);
        }
    };
}

#endif
