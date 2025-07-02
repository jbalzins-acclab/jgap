#ifndef EAMKERNEL_HPP
#define EAMKERNEL_HPP

#include "core/kernels/Kernel.hpp"
#include "data/kernels/EamKernelIndex.hpp"



namespace jgap {

    class EamKernelSE : public Kernel<double, EamKernelIndex> {
    public:

        EamKernelSE(const double energyScale, const double lengthScale) : _lengthScale(lengthScale) {
            _inverse2ThetaSq = 1.0 / (2.0 * _lengthScale * _lengthScale);
            _energyScaleSquared = energyScale * energyScale;
        }

        ~EamKernelSE() override = default;

        double covariance(
            const AtomicStructure &structure,
            const EamKernelIndex &indexes,
            const double &sparseDensity) override;

        vector<Vector3> derivatives(
            const AtomicStructure &structure,
            const EamKernelIndex &indexes,
            const double &sparseDensity) override;

        [[nodiscard]] double covariance(const double &density1, const double &density2) override;

    private:

        [[nodiscard]] double derivative(const double &density1, const double &density2);

        double _energyScaleSquared;
        double _lengthScale;
        double _inverse2ThetaSq;
    };

}

#endif //EAMKERNEL_HPP
