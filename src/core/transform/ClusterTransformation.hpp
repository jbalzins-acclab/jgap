#ifndef JGAP_TRANSFORMER_HPP
#define JGAP_TRANSFORMER_HPP

#include <vector>
#include <map>

#include "../atomic/Descriptor.hpp"
#include "core/atomic/geometry/Cluster.hpp"
#include "core/atomic/neighbours/NeighbourList.hpp"
#include "../potentials/Cutoffs.hpp"
#include "io/Serializable.hpp"
#include "utils/ValuePtr.hpp"

namespace jgap {

    template<size_t NDim, size_t NClusterSize>
    class ClusterTransformation {
    public:
        static constexpr size_t Dim = NDim;
        static constexpr size_t ClusterSize = NClusterSize;

        virtual ~ClusterTransformation() = default;

        virtual Descriptor<Dim> evaluate(const Cluster<ClusterSize>& cluster) const = 0;
        virtual NBodyDescriptor<Dim, ClusterSize> evaluateAndDifferentiate(const Cluster<ClusterSize>& cluster)
            const = 0;

        virtual Cutoffs getCutoffs() const = 0;

        virtual Real symmetryFactor() const {
            return 1.0;
        }

        virtual std::unique_ptr<ClusterTransformation> clone() const = 0;
    };

    static_assert(Cloneable<ClusterTransformation<1, 2>>);
}

#endif