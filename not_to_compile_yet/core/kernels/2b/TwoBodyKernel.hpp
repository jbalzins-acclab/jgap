#ifndef JGAP_TWOBODYKERNEL_HPP
#define JGAP_TWOBODYKERNEL_HPP

#include "../Kernel.hpp"
#include "core/cutoff/CutoffFunction.hpp"
#include "../../../data/atomic/AtomicStructure.hpp"
#include "../../../data/atomic/PredictionData.hpp"
#include "data/descriptors/kernels/TwoBodyIndex.hpp"
namespace jgap {

    class TwoBodyKernel : public Kernel<SpeciesPair, TwoBodyIndex, TwoBodyDescriptorData> {
    public:
        Covariance covariance(const AtomicStructure &structure, const TwoBodyIndex &index) override;
        double value(const TwoBodyDescriptorData &r) override;
    protected:
        virtual double valueInternal(const double &r) const = 0;
        virtual double derivativeInternal(const double &changingR) const = 0;
    };
}

#endif