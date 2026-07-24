#ifndef JGAP_TWOBODYTRANSFORMATION_HPP
#define JGAP_TWOBODYTRANSFORMATION_HPP

#include "jgap/core/ValuePtr.hpp"
#include "jgap/core/atomic/descriptor/Descriptor.hpp"
#include "jgap/core/atomic/descriptor/TwoBodyDescriptor.hpp"
#include "jgap/core/atomic/geometry/Cluster2.hpp"
#include "jgap/core/potentials/Cutoffs.hpp"

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
