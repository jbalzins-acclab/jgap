#ifndef JGAP_TRANSFORMER_HPP
#define JGAP_TRANSFORMER_HPP

#include <vector>
#include <map>

#include "../atomic/Descriptor.hpp"
#include "core/atomic/geometry/Separations.hpp"
#include "core/atomic/neighbours/NeighbourList.hpp"
#include "data/Cutoffs.hpp"
#include "io/Serializable.hpp"

namespace jgap {

    template<size_t Dim, size_t DependencyAtoms = GradientData::UNKNOWN_DEPENDENCIES>
    class Transformer {// : public Serializable {
    public:
        using TValueOnly = Descriptor<Dim, DependencyAtoms, CalculationType::ValueOnly>;
        using TWithGradients = Descriptor<Dim, DependencyAtoms, CalculationType::WithGradients>;

        virtual ~Transformer() = default;

        virtual Real getCutoff() const = 0;

        virtual std::map<SpeciesSet, std::vector<TValueOnly>> transform(
            const NeighbourList& neighbour_list) const = 0;
        virtual std::map<SpeciesSet, std::vector<TWithGradients>> transformWithGradients(
            const NeighbourList& neighbour_list) const = 0;

        Cutoffs getCutoffs() {
            return Cutoffs{{{DependencyAtoms, getCutoff()}}};
        }
    };
}

#endif