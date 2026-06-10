#ifndef JGAP_TRANSFORMATIONAGGREGATOR_HPP
#define JGAP_TRANSFORMATIONAGGREGATOR_HPP

#include "core/atomic/neighbours/NeighbourList.hpp"
#include "core/transform/ClusterTransformation.hpp"
#include "core/atomic/species/Species.hpp"
#include "core/atomic/species/SpeciesSet.hpp"
#include "io/log/CurrentLogger.hpp"
#include <map>
#include <vector>
#include <memory>
#include <array>

namespace jgap {

    template<size_t Dim>
    class TransformationAggregator {
    public:
        virtual ~TransformationAggregator() = default;
        virtual std::map<size_t, ManyBodyDescriptor<Dim>> aggregate(const NeighbourList& nl) const = 0;
        virtual Cutoffs getCutoffs() const = 0;
    };
}

#endif