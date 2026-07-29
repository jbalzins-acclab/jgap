#ifndef JGAP_KERNEL_HPP
#define JGAP_KERNEL_HPP

#include <array>
#include <memory>
#include <optional>
#include <string>
#include <vector>
#include "jgap/core/ValuePtr.hpp"
#include "jgap/core/atomic/descriptor/Descriptor.hpp"

namespace jgap {

    template<size_t NDim>
    class Kernel {
    public:
        constexpr static size_t Dim = NDim;

        struct KernelValueAndGradient {
            Real value;
            std::array<Real, Dim> gradient;
        };

        virtual ~Kernel() = default;

        /// @note In derived classes, overriding value by calling Kernel::value(q1, q2)
        /// (e.g. `Real value(const Descriptor<Dim>& q1, const Descriptor<Dim>& q2) const override { return Kernel::value(q1, q2); }`)
        /// forces devirtualization of valueAndGradient for compiler optimizations.
        virtual Real value(const Descriptor<Dim>& q1, const Descriptor<Dim>& q2) const {
            return valueAndGradient(q1, q2).value;
        }

        virtual KernelValueAndGradient valueAndGradient(
            const Descriptor<Dim>& sparse_point, const Descriptor<Dim>& q
        ) const = 0;

        virtual Kernel<NDim>* clone() const = 0;
    };

    static_assert(Cloneable<Kernel<1>>);
    static_assert(Cloneable<Kernel<2>>);

    template<typename T, size_t Dim>
    concept CKernelOfDim = std::derived_from<T, Kernel<Dim>>;
}

#endif
