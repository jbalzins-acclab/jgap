#ifndef JGAP_THREEBODYTRANSFORMATION_HPP
#define JGAP_THREEBODYTRANSFORMATION_HPP

#include "jgap/core/ValuePtr.hpp"
#include "jgap/core/atomic/descriptor/ThreeBodyDescriptor.hpp"
#include "jgap/core/atomic/geometry/Cluster3.hpp"
#include "jgap/core/potentials/Cutoffs.hpp"

namespace jgap {
    template<size_t Dim>
    class ThreeBodyTransformation {
    public:
        virtual ~ThreeBodyTransformation() = default;

        /// @note In derived classes, overriding evaluate by calling ThreeBodyTransformation::evaluate(cluster)
        /// (e.g. `Descriptor<Dim> evaluate(const Cluster3& cluster) const override { return ThreeBodyTransformation::evaluate(cluster); }`)
        /// forces devirtualization of evaluateAndDifferentiate for compiler optimizations.
        virtual Descriptor<Dim> evaluate(const Cluster3& cluster) const {
            return evaluateAndDifferentiate(cluster).value;
        }

        virtual ThreeBodyDescriptor<Dim> evaluateAndDifferentiate(const Cluster3& cluster) const = 0;
        virtual Cutoffs getCutoffs() const = 0;
        virtual bool isRotationallyInvariant() const = 0;
        virtual bool isSwapInvariant(size_t idx1, size_t idx2) const = 0;
        virtual ThreeBodyTransformation<Dim>* clone() const = 0;
    };

    static_assert(Cloneable<ThreeBodyTransformation<1>>);
    static_assert(Cloneable<ThreeBodyTransformation<2>>);
}
#endif
