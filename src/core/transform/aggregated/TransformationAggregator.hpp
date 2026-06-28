#ifndef JGAP_TRANSFORMATIONAGGREGATOR_HPP
#define JGAP_TRANSFORMATIONAGGREGATOR_HPP

#include "core/atomic/neighbours/NeighbourList.hpp"
#include "core/transform/NBodyTransformation.hpp"
#include "core/atomic/species/Species.hpp"
#include "core/atomic/species/SpeciesSet.hpp"
#include "io/log/CurrentLogger.hpp"
#include <map>
#include <vector>
#include <memory>
#include <array>

#include "core/tabulation/TabulationData.hpp"

namespace jgap {

    template<size_t Dim>
    class TransformationAggregator {
    public:
        virtual ~TransformationAggregator() = default;
        virtual ManyBodyDescriptor<Dim> aggregate(const NeighbourList& nl) const = 0;
        virtual Cutoffs getCutoffs() const = 0;
        virtual TransformationAggregator* clone() const = 0;
        virtual Species getCentralSpecies() const = 0;
        virtual std::set<Species> getAllSpecies() const = 0;
        virtual void tabulateNewManyBodyGrid(TabulationData& tables) const = 0;
    };

    static_assert(Cloneable<TransformationAggregator<1>>);
    static_assert(Cloneable<TransformationAggregator<2>>);
    static_assert(Cloneable<TransformationAggregator<3>>);
}

#endif