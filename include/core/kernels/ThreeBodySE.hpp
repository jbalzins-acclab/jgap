#ifndef THREEBODYSE_HPP
#define THREEBODYSE_HPP

#include "../cutoff/CutoffFunction.hpp"
#include "core/kernels/Kernel.hpp"
#include "data/kernels/ThreeBodyKernelIndex.hpp"

namespace jgap {
    class ThreeBodySE : public Kernel<Vector3, ThreeBodyKernelIndex> {
    public:
        ThreeBodySE(const shared_ptr<CutoffFunction> &cutoffFunction,
            const double energyScale, const double lengthScale)
            : _cutoffFunction(cutoffFunction), _lengthScale(lengthScale) {
            _inverse2SigmaSq = 1.0 / (2.0 * _lengthScale * _lengthScale);
            _energyScaleSquared = energyScale * energyScale;
        }

        ~ThreeBodySE() override = default;

        double covariance(const AtomicStructure &structure,
                          const ThreeBodyKernelIndex &indexes,
                          const Vector3 &descriptorInvariantDistances) override;
        vector<Vector3> derivatives(const AtomicStructure &structure,
                                    const ThreeBodyKernelIndex &indexes,
                                    const Vector3 &descriptorInvariantDistances) override;

        double covariance(const Vector3 &t1, const Vector3 &t2) override;

    private:
        shared_ptr<CutoffFunction> _cutoffFunction;
        double _energyScaleSquared;
        double _lengthScale;
        double _inverse2SigmaSq;

        [[nodiscard]] Vector3 gradient(const Vector3 &changingTriplet, const Vector3 &constTriplet);
    };
}

#endif //THREEBODYSE_HPP
