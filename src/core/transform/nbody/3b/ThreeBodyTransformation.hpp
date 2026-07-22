#ifndef JGAP_THREEBODYTRANSFORMATION_HPP
#define JGAP_THREEBODYTRANSFORMATION_HPP

#include "core/ValuePtr.hpp"
#include "core/atomic/descriptor/ThreeBodyDescriptor.hpp"
#include "core/atomic/geometry/Cluster3.hpp"
#include "core/potentials/Cutoffs.hpp"

namespace jgap {
    template<size_t Dim>
    class ThreeBodyTransformation {
    public:
        virtual ~ThreeBodyTransformation() = default;
        virtual Descriptor<Dim> evaluate(const Cluster3& cluster) const = 0;
        virtual ThreeBodyDescriptor<Dim> evaluateAndDifferentiate(const Cluster3& cluster) const = 0;
        virtual Cutoffs getCutoffs() const = 0;
        virtual bool isSwapInvariant(size_t idx1, size_t idx2) const = 0;
        virtual ThreeBodyTransformation<Dim>* clone() const = 0;
    };

    static_assert(Cloneable<ThreeBodyTransformation<1>>);
    static_assert(Cloneable<ThreeBodyTransformation<2>>);
}
#endif
