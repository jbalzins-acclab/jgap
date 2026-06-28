#ifndef JGAP_NBODYTRANSFORMATION_HPP
#define JGAP_NBODYTRANSFORMATION_HPP

#include <vector>
#include <map>

#include "../atomic/Descriptor.hpp"
#include "core/atomic/geometry/Cluster.hpp"
#include "core/atomic/neighbours/NeighbourList.hpp"
#include "../potentials/Cutoffs.hpp"
#include "../ValuePtr.hpp"
#include "core/atomic/iteration/ClusterFinder.hpp"

namespace jgap {

    template<size_t NDim, size_t NClusterSize>
    class NBodyTransformation {
    public:
        static constexpr size_t Dim = NDim;
        static constexpr size_t ClusterSize = NClusterSize;

        virtual ~NBodyTransformation() = default;

        virtual Descriptor<Dim> evaluate(const Cluster<ClusterSize>& cluster) const = 0;
        virtual NBodyDescriptor<Dim, ClusterSize> evaluateAndDifferentiate(const Cluster<ClusterSize>& cl) const = 0;

        virtual Cutoffs getCutoffs() const = 0;


        virtual bool isSymmetricWrtNodeIndices() const = 0;

        virtual NBodyTransformation* clone() const = 0;

    };

    static_assert(Cloneable<NBodyTransformation<1, 2>>);
}

#endif