#ifndef JGAP_EAMKERNEL_HPP
#define JGAP_EAMKERNEL_HPP

#include "core/descriptors/Kernel.hpp"
#include "data/descriptors/kernels/EamKernelIndex.hpp"

namespace jgap {
    class EamKernel : public Kernel<Species, EamKernelIndexPerSpecies, double> {
    public:
        Covariance covariance(const AtomicStructure &structure, const EamKernelIndexPerSpecies &indexes) override;
        double value(const double &density) override;

        // TODO: optional?
        virtual double getDensity() = 0;

    protected:
        virtual double valueInternal(const double &density) const = 0;
        virtual double derivativeInternal(const double &density) const = 0;
    };
}

#endif