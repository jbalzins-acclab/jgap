#ifndef TWOBODYSE_HPP
#define TWOBODYSE_HPP

#include "core/kernels/Kernel.hpp"
#include "../cutoff/CutoffFunction.hpp"
#include "data/kernels/TwoBodyKernelIndex.hpp"

namespace jgap {
    class TwoBodySE : public Kernel<double, TwoBodyKernelIndex> {
    public:
        TwoBodySE(const shared_ptr<CutoffFunction> &cutoffFunction, const double energyScale, const double lengthScale)
            : _cutoffFunction(cutoffFunction), _lengthScale(lengthScale) {
            _inverse2ThetaSq = 1.0 / (2.0 * _lengthScale * _lengthScale);
            _energyScaleSquared = energyScale * energyScale;
        }
        ~TwoBodySE() override = default;

        double covariance(const AtomicStructure &structure,
                          const TwoBodyKernelIndex &indexes,
                          const double &rSparse) override;
        vector<Vector3> derivatives(const AtomicStructure &structure,
                                    const TwoBodyKernelIndex &indexes,
                                    const double &rSparse) override;

        double covariance(const double &r1, const double &r2) override;

    private:
        shared_ptr<CutoffFunction> _cutoffFunction;
        double _energyScaleSquared;
        double _lengthScale;
        double _inverse2ThetaSq;

        // TODO: is this the best way to organize it?
        double derivative(const double &changingR, const double &constR);
    };
}

#endif //TWOBODYSE_HPP
