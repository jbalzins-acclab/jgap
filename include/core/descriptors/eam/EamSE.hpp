#ifndef EAMKERNEL_HPP
#define EAMKERNEL_HPP

#include "EamKernel.hpp"
#include "../Kernel.hpp"
#include "data/descriptors/kernels/EamKernelIndex.hpp"
#include "io/parse/ParserRegistry.hpp"

namespace jgap {

    class EamSE : public EamKernel {
    public:
        static constexpr string TYPE = "squared_exp";
        EamSE(Species species, double energyScale, double lengthScale, double density, optional<double> coeff = {});

        EamSE(const nlohmann::json &params);
        string getType() override { return TYPE; }
        nlohmann::json serialize() override;

        double crossCovariance(const shared_ptr<IKernel>& other) override;

        Species getFilter() override { return _idSpecies; }
        pair<double, double> getDensityRange() override { return {_density, _density}; }

    private:
        // raw params
        Species _idSpecies;
        double _energyScale;
        double _lengthScale;
        double _density;

        // optimized for calculation
        double _totalPrefactor;
        double _inverseThetaSq;

        double valueInternal(const double &density) const override;
        double derivativeInternal(const double &density) const override;
    };

    REGISTER_PARSER(EamKernel, EamSE)
}

#endif
