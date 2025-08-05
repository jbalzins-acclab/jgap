#ifndef THREEBODYSE_HPP
#define THREEBODYSE_HPP

#include "core/cutoff/CutoffFunction.hpp"
#include "core/descriptors/kernels/Kernel.hpp"
#include "data/descriptors/kernels/ThreeBodyKernelIndex.hpp"
#include "io/parse/ParserRegistry.hpp"

#include <nlohmann/json.hpp>

#include "io/parse/ParserRegistry.hpp"

namespace jgap {

    using ThreeBodyKernel = Kernel<ThreeBodyDescriptorData, ThreeBodyKernelIndex>;

    class ThreeBodySE : public ThreeBodyKernel {
    public:
        explicit ThreeBodySE(double energyScale, double lengthScale);

        explicit ThreeBodySE(const nlohmann::json& params);
        string getType() override { return "squared_exp"; }
        nlohmann::json serialize() override;

        Covariance covariance(const AtomicStructure &structure,
                              const ThreeBodyKernelIndex &indexes,
                              const ThreeBodyDescriptorData &descriptorInvariantDistances) override;

        double covariance(const ThreeBodyDescriptorData &t1, const ThreeBodyDescriptorData &t2) override;

    private:
        double _energyScaleSquared;
        double _lengthScale;
        double _inverse2ThetaSq;
        double _inverseThetaSq;

        [[nodiscard]] double covarianceNoCutoffs(const Vector3 &t1, const Vector3 &t2) const;
        [[nodiscard]] Vector3 gradient(const Vector3 &changingQ, const Vector3 &constQ) const;
    };

    REGISTER_PARSER("squared_exp", ThreeBodyKernel, ThreeBodySE);
}

#endif //THREEBODYSE_HPP
