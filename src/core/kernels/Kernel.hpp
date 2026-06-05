#ifndef JGAP_KERNEL_HPP
#define JGAP_KERNEL_HPP

#include "../atomic/geometry/Vector3.hpp"
#include <optional>
#include <string>
#include <memory>
#include <array>
#include <vector>

#include "../atomic/Descriptor.hpp"
#include "core/atomic/energy/AtomicQuantity.hpp"
#include "core/atomic/geometry/Cluster.hpp"

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

        virtual Real value(const std::array<Real, Dim>& q1, const std::array<Real, Dim>& q2) const = 0;

        virtual KernelValueAndGradient valueAndGradient(const std::array<Real, Dim> &sparse_point,
                                                        const std::array<Real, Dim> &q) const = 0;

    };

    template<typename T>
    concept CKernel = std::derived_from<T, Kernel<T::Dim>>;
}

#endif
