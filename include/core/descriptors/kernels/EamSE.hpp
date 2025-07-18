#ifndef EAMKERNEL_HPP
#define EAMKERNEL_HPP

#include "core/descriptors/kernels/Kernel.hpp"
#include "data/kernels/EamKernelIndex.hpp"
#include "io/parse/ParserRegistry.hpp"

namespace jgap {

    using EamKernel = Kernel<double, EamKernelIndexPerSpecies>;

    class EamSE : public EamKernel {
    public:

        explicit EamSE(const nlohmann::json &params);
        explicit EamSE(double energyScale, double lengthScale);

        Covariance covariance(const AtomicStructure &structure,
                              const EamKernelIndexPerSpecies &indexes,
                              const double &sparseDensity) override;

        double covariance(const double &density1, const double &density2) override;

        string getType() override { return "squared_exp"; }
        nlohmann::json serialize() override;

    private:
        double _energyScaleSquared;
        double _lengthScale;
        double _inverse2ThetaSq;

        [[nodiscard]] double covarianceNoCutoffs(const double &density1, const double &density2) const;
        [[nodiscard]] double derivative(const double &changingDensity, const double &constantDensity) const;
    };

    REGISTER_PARSER("squared_exp", EamKernel, EamSE)
}

#endif //EAMKERNEL_HPP
