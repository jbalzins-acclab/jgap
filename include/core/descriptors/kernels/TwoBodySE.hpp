#ifndef TWOBODYSE_HPP
#define TWOBODYSE_HPP

#include "core/descriptors/kernels/Kernel.hpp"
#include "core/cutoff/CutoffFunction.hpp"
#include "data/kernels/TwoBodyKernelIndex.hpp"
#include "io/parse/ParserRegistry.hpp"

namespace jgap {

    using TwoBodyKernel = Kernel<double, TwoBodyKernelIndex>;

    class TwoBodySE : public TwoBodyKernel {
    public:
        explicit TwoBodySE(const nlohmann::json &params);
        explicit TwoBodySE(const shared_ptr<CutoffFunction> &cutoffFunction, double energyScale, double lengthScale);

        ~TwoBodySE() override = default;

        Covariance covariance(const AtomicStructure &structure,
                              const TwoBodyKernelIndex &indexes,
                              const double &rSparse) override;

        double covariance(const double &r1, const double &r2) override;

        string getType() override { return "squared_exp"; }
        nlohmann::json serialize() override;

    private:
        shared_ptr<CutoffFunction> _cutoffFunction;
        double _energyScaleSquared;
        double _lengthScale;
        double _inverse2ThetaSq;

        double covarianceNoCutoffs(const double &r1, const double &r2) const;
        double derivativeNoCutoffs(const double &changingR, const double &constR) const;
    };

    REGISTER_PARSER("squared_exp", TwoBodyKernel, TwoBodySE)
}

#endif //TWOBODYSE_HPP
