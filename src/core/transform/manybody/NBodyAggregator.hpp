#ifndef JGAP_NBODYAGGREGATOR_HPP
#define JGAP_NBODYAGGREGATOR_HPP

#include <array>
#include <map>
#include <memory>
#include <vector>
#include "core/ValuePtr.hpp"
#include "core/atomic/descriptor/ManyBodyDescriptors.hpp"
#include "core/atomic/neighbours/NeighbourLists.hpp"
#include "core/atomic/species/Species.hpp"

#include "core/tabulation/TabulationData.hpp"

namespace jgap {

    template<size_t Dim>
    class NBodyAggregator {
    public:
        virtual ~NBodyAggregator() = default;
        virtual ManyBodyDescriptors<Dim> aggregate(const NeighbourLists& nl) const = 0;
        virtual Cutoffs getCutoffs() const = 0;
        virtual NBodyAggregator* clone() const = 0;
        virtual Species getCentralSpecies() const = 0;
        virtual std::set<Species> getAllSpecies() const = 0;
        virtual void tabulateNewManyBodyGrid(TabulationData& tables) const = 0;
    };

    static_assert(Cloneable<NBodyAggregator<1>>);
    static_assert(Cloneable<NBodyAggregator<2>>);
    static_assert(Cloneable<NBodyAggregator<3>>);
    static_assert(Cloneable<NBodyAggregator<4>>);
}

#endif
