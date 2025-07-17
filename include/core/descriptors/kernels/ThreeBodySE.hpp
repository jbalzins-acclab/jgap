#ifndef THREEBODYSE_HPP
#define THREEBODYSE_HPP

#include "core/cutoff/CutoffFunction.hpp"
#include "core/descriptors/kernels/Kernel.hpp"
#include "data/kernels/ThreeBodyKernelIndex.hpp"
#include "io/parse/ParserRegistry.hpp"

#include <nlohmann/json.hpp>

#include "io/parse/ParserRegistry.hpp"

namespace jgap {

    using ThreeBodyKernel = Kernel<Vector3, ThreeBodyKernelIndex>;

    class ThreeBodySE : public ThreeBodyKernel {
    public:
        explicit ThreeBodySE(const nlohmann::json& params);
        explicit ThreeBodySE(const shared_ptr<CutoffFunction> &cutoffFunction, double energyScale, double lengthScale);

        ~ThreeBodySE() override = default;

        double covariance(const AtomicStructure &structure,
                          const ThreeBodyKernelIndex &indexes,
                          const Vector3 &descriptorInvariantDistances) override;
        vector<Vector3> derivatives(const AtomicStructure &structure,
                                    const ThreeBodyKernelIndex &indexes,
                                    const Vector3 &descriptorInvariantDistances) override;

        double covariance(const Vector3 &t1, const Vector3 &t2) override;

        string getType() override { return "squared_exp"; }
        nlohmann::json serialize() override;

    private:
        double _energyScaleSquared;
        double _lengthScale;
        double _inverse2ThetaSq;

        shared_ptr<CutoffFunction> _cutoffFunction;

        [[nodiscard]] double invariantTripletToCutoff(const Vector3 &t) const;

        [[nodiscard]] double covarianceNoCutoffs(const Vector3 &t1, const Vector3 &t2) const;
        [[nodiscard]] Vector3 gradient(const Vector3 &changingTriplet, const Vector3 &constTriplet) const;
    };

    REGISTER_PARSER("squared_exp", ThreeBodyKernel, ThreeBodySE);
}

#endif //THREEBODYSE_HPP
