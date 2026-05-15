#ifndef JGAP_KNOWNDEPENDENCIESTRANSFORMER_HPP
#define JGAP_KNOWNDEPENDENCIESTRANSFORMER_HPP

#include "Transformer.hpp"
#include "core/atomic/Descriptor.hpp"

namespace jgap {
    template<size_t Dim, size_t DependencyAtoms>
    requires (DependencyAtoms > 0)
    class KnownDependenciesTransformer : public Transformer<Dim, DependencyAtoms> {
    public:
        using TValueOnly = Descriptor<Dim, DependencyAtoms, CalculationType::ValueOnly>;
        using TWithGradients = Descriptor<Dim, DependencyAtoms, CalculationType::WithGradients>;

        using TValue = std::array<Real, Dim>;
        using TDerivatives = std::array<std::array<Real, Dim>, Separations<DependencyAtoms>::N_SEPARATIONS>;

        struct ValueAndDerivatives {
            TValue value;
            TDerivatives derivatives;
        };

        std::map<SpeciesSet, std::vector<TValueOnly>> transform(
            const NeighbourList& neighbour_list) const override;

        std::map<SpeciesSet, std::vector<TWithGradients>> transformWithGradients(
            const NeighbourList& neighbour_list) const override;

        template<CalculationType Type>
        Descriptor<Dim, DependencyAtoms, Type> transform(const Separations<DependencyAtoms>& separations) const;

        virtual TValue value(const Separations<DependencyAtoms>& separations) const = 0;

        // d\vec{q}/d|r_ij|
        virtual ValueAndDerivatives valueAndDerivatives(const Separations<DependencyAtoms>& separations) const = 0;
    };

    template<size_t Dim, size_t DependencyAtoms> requires (DependencyAtoms > 0)
    std::map<SpeciesSet, std::vector<typename KnownDependenciesTransformer<Dim, DependencyAtoms>::TValueOnly>>
    KnownDependenciesTransformer<Dim, DependencyAtoms>::transform(const NeighbourList &neighbour_list) const {

        double cutoff = this->getCutoff();

        std::map<SpeciesSet, std::vector<TValueOnly>> result;

        for (const auto& species_set: neighbour_list.getSpeciesSets<DependencyAtoms>()) {
            result[species_set] = {};
            for (const auto& separations: neighbour_list.findAllSeparations<DependencyAtoms>(species_set, cutoff)) {
                result[species_set].push_back(this->template transform<CalculationType::ValueOnly>(separations));
            }
        }

        return result;
    }

    template<size_t Dim, size_t DependencyAtoms> requires (DependencyAtoms > 0)
    std::map<SpeciesSet, std::vector<typename KnownDependenciesTransformer<Dim, DependencyAtoms>::TWithGradients>>
    KnownDependenciesTransformer<Dim, DependencyAtoms>::transformWithGradients(
        const NeighbourList &neighbour_list) const {

        double cutoff = this->getCutoff();

        std::map<SpeciesSet, std::vector<TWithGradients>> result;

        for (const auto& species_set: neighbour_list.getSpeciesSets<DependencyAtoms>()) {
            result[species_set] = {};
            for (const auto& separations: neighbour_list.findAllSeparations<DependencyAtoms>(species_set, cutoff)) {
                result[species_set].push_back(this->template transform<CalculationType::WithGradients>(separations));
            }
        }

        return result;
    }

    template<size_t Dim, size_t DependencyAtoms> requires (DependencyAtoms > 0)
    template<CalculationType Type>
    Descriptor<Dim, DependencyAtoms, Type> KnownDependenciesTransformer<Dim, DependencyAtoms>::transform(
        const Separations<DependencyAtoms> &separations) const {

        if constexpr (Type == CalculationType::ValueOnly) {

            return Descriptor<Dim, DependencyAtoms>{ value(separations) };

        } else {

            auto [value, derivatives_wrt_r_norms] = valueAndDerivatives(separations);

            Descriptor<Dim, DependencyAtoms, Type> q_and_derivatives;

            for (size_t dim = 0; dim < Dim; dim++) {
                q_and_derivatives.virials[dim] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
                for (size_t i = 0; i < DependencyAtoms; i++) {
                    q_and_derivatives.gradients[dim][i].wrt_atom_index = separations.atom_indexes[i];
                }
            }

            //auto  = derivatives(separations);
            for (size_t i = 0; i < DependencyAtoms; i++) {
                for (size_t j = i + 1; j < DependencyAtoms; j++) {
                    const auto sep_idx = flattened_index(i, j);
                    const auto& sep = separations.between(i, j);

                    for (size_t dim = 0; dim < Dim; dim++) {
                        const auto deriv = derivatives_wrt_r_norms[sep_idx][dim];

                        q_and_derivatives.gradients[dim][j].value += sep.direction * deriv;
                        q_and_derivatives.gradients[dim][i].value -= sep.direction * deriv;

                        //TODO: JUST STORE SEPARATION VIRIALS,
                        q_and_derivatives.virials[dim] += sep.virials * deriv;
                    }
                }
            }

            return q_and_derivatives;
        }
    }
}

#endif