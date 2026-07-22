#ifndef JGAP_KERNEL_HPP
#define JGAP_KERNEL_HPP

#include <array>
#include <memory>
#include <optional>
#include <string>
#include <vector>
#include "core/atomic/descriptor/Descriptor.hpp"

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

        virtual Real value(const Descriptor<Dim> &q1, const Descriptor<Dim> &q2) const = 0;

        virtual KernelValueAndGradient valueAndGradient(const Descriptor<Dim> &sparse_point,
                                                        const Descriptor<Dim> &q) const = 0;
    };

    template<typename T>
    concept CKernel = std::derived_from<T, Kernel<T::Dim>>;

    template<typename T, size_t Dim>
    concept CKernelOfDim = std::derived_from<T, Kernel<Dim>>;
}

#endif
