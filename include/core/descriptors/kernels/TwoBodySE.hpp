#ifndef TWOBODYSE_HPP
#define TWOBODYSE_HPP

#include "core/descriptors/kernels/Kernel.hpp"
#include "core/cutoff/CutoffFunction.hpp"
#include "data/descriptors/kernels/TwoBodyKernelIndex.hpp"
#include "io/parse/ParserRegistry.hpp"

namespace jgap {

    using TwoBodyKernel = Kernel<TwoBodyDescriptorData, TwoBodyKernelIndex>;

    class TwoBodySE : public TwoBodyKernel {
    public:
        explicit TwoBodySE(const nlohmann::json &params);
        string getType() override { return "squared_exp"; }
        nlohmann::json serialize() override;

        Covariance covariance(const AtomicStructure &structure,
                              const TwoBodyKernelIndex &indexes,
                              const TwoBodyDescriptorData &rSparse) override;

        double covariance(const TwoBodyDescriptorData &r1, const TwoBodyDescriptorData &r2) override;

    private:
        double _energyScaleSquared;
        double _lengthScale;
        double _inverse2ThetaSq;
        double _inverseThetaSq;

        double covarianceNoCutoffs(const double &r1, const double &r2) const;
        double derivativeNoCutoffs(const double &changingR, const double &constR) const;
    };

    REGISTER_PARSER("squared_exp", TwoBodyKernel, TwoBodySE)
}

#endif //TWOBODYSE_HPP
