#ifndef JGAP_CUMULATIVETRANSFORMER_HPP
#define JGAP_CUMULATIVETRANSFORMER_HPP

#include "Transformer.hpp"
#include <map>
#include <vector>

namespace jgap {

    template<size_t Dim, size_t DependencyAtoms/*PerContributor*/>
    class CumulativePerAtomTransformer : public Transformer<Dim> {
    public:
        using TValueOnly = Descriptor<Dim, GradientData::UNKNOWN_DEPENDENCIES, CalculationType::ValueOnly>;
        using TWithGradients = Descriptor<Dim, GradientData::UNKNOWN_DEPENDENCIES, CalculationType::WithGradients>;

        std::map<SpeciesSet, std::vector<TValueOnly>> transform(
            const NeighbourList& neighbour_list) const override {

            std::map<SpeciesSet, std::vector<TValueOnly>> result_map;
            std::vector<TValueOnly> result(neighbour_list.nAtoms());

            double cutoff = this->getCutoff();

            for (const auto& species_set: neighbour_list.getSpeciesSets<DependencyAtoms>()) {
                for (const auto& separations: neighbour_list.findAllSeparations<DependencyAtoms>(species_set, cutoff)) {
                    contribution<CalculationType::ValueOnly>(separations, result[separations.atom_indexes[0]]);
                }
            }

            // Assuming cumulative transformers map the central atom's species to its descriptors
            for (size_t i = 0; i < neighbour_list.nAtoms(); i++) {
                SpeciesSet central_species{neighbour_list.speciesOf(i)};
                result_map[central_species].push_back(result[i]);
            }

            return result_map;
        }

        std::map<SpeciesSet, std::vector<TWithGradients>> transformWithGradients(
            const NeighbourList& neighbour_list) const override {

            std::map<SpeciesSet, std::vector<TWithGradients>> result_map;
            std::vector<TWithGradients> result(neighbour_list.nAtoms(), TWithGradients(neighbour_list.nAtoms()));

            double cutoff = this->getCutoff();

            for (const auto& species_set: neighbour_list.getSpeciesSets<DependencyAtoms>()) {
                for (const auto& separations: neighbour_list.findAllSeparations<DependencyAtoms>(species_set, cutoff)) {
                    contribution<CalculationType::WithGradients>(separations, result[separations.atom_indexes[0]]);
                }
            }

            for (size_t i = 0; i < neighbour_list.nAtoms(); i++) {
                SpeciesSet central_species{neighbour_list.speciesOf(i)};
                result_map[central_species].push_back(result[i]);
            }

            return result_map;
        }

        template<CalculationType Type>
        void contribution(const Separations<DependencyAtoms>& separations, Descriptor<Dim>& contributed_to) {

            std::array<Real, Dim> q = value(separations);

            for (int i = 0; i < Dim; i++) {
                contributed_to.value[i] += q[i];
            }

            if constexpr (Type == CalculationType::ValueOnly) {
                return;
            }

            auto dq_drij = derivatives(separations);
            for (size_t i = 0; i < DependencyAtoms; i++) {
                for (size_t j = i + 1; j < DependencyAtoms; j++) {
                    const auto sep_idx = flattened_index(i, j);
                    const auto& sep = separations.between(i, j);
                    for (int dim = 0; dim < Dim; dim++) {
                        contributed_to.virials[dim] += sep.virials * dq_drij[sep_idx][dim];
                        contributed_to.gradients[dim][sep.atom_indexes[j]] += sep.direction * dq_drij[sep_idx][dim];
                        contributed_to.gradients[dim][sep.atom_indexes[i]] -= sep.direction * dq_drij[sep_idx][dim];
                    }
                }
            }
        }

        virtual std::array<Real, Dim> value(const Separations<DependencyAtoms>& separations) const = 0;
        virtual std::array<std::array<Real, Dim>, Separations<DependencyAtoms>::N_SEPARATIONS> derivatives(
            const Separations<DependencyAtoms>& separations) const = 0;
    };
}

#endif
