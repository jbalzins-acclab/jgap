#ifndef JGAP_THREEBODYKERNEL_HPP
#define JGAP_THREEBODYKERNEL_HPP

#include "../Kernel.hpp"
#include "data/descriptors/kernels/ThreeBodyIndex.hpp"

namespace jgap {
    class ThreeBodyKernel : public Kernel<SpeciesTriplet, ThreeBodyIndex, ThreeBodyDescriptorData> {
    public:
        Covariance covariance(const AtomicStructure& structure, const ThreeBodyIndex& index) override;
        double value(const ThreeBodyDescriptorData& t) override;
    protected:
        virtual double valueInternal(const Vector3& q) const = 0;
        virtual Vector3 gradientInternal(const Vector3& q) const = 0;
    };
}

#endif