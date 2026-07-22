#ifndef JGAP_TWOBODYTRANSFORMATION_HPP
#define JGAP_TWOBODYTRANSFORMATION_HPP

#include "core/ValuePtr.hpp"
#include "core/atomic/descriptor/Descriptor.hpp"
#include "core/atomic/descriptor/TwoBodyDescriptor.hpp"
#include "core/atomic/geometry/Cluster2.hpp"
#include "core/potentials/Cutoffs.hpp"

namespace jgap {
    template<size_t Dim>
    class TwoBodyTransformation {
    public:
        virtual ~TwoBodyTransformation() = default;
        virtual Descriptor<Dim> evaluate(const Cluster2& pair) const = 0;
        virtual TwoBodyDescriptor<Dim> evaluateAndDifferentiate(const Cluster2& pair) const = 0;
        virtual Cutoffs getCutoffs() const = 0;
        virtual TwoBodyTransformation<Dim>* clone() const = 0;
    };

    static_assert(Cloneable<TwoBodyTransformation<1>>);
    static_assert(Cloneable<TwoBodyTransformation<2>>);
}
#endif
