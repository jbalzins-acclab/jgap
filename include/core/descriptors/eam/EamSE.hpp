#ifndef EAMKERNEL_HPP
#define EAMKERNEL_HPP

#include "../Kernel.hpp"
#include "data/descriptors/kernels/EamKernelIndex.hpp"
#include "io/parse/ParserRegistry.hpp"

namespace jgap {

    using EamKernel = Kernel<Species, EamKernelIndexPerSpecies, double>;

    class EamSE : public EamKernel {
    public:
        EamSE(Species species, double energyScale, double lengthScale, double density, optional<double> coeff = {});

        EamSE(const nlohmann::json &params);
        string getType() override { return "squared_exp"; }
        nlohmann::json serialize() override;

        Covariance covariance(const AtomicStructure &structure, const EamKernelIndexPerSpecies &indexes) override;
        double crossCovariance(const shared_ptr<IKernel>& other) override;

        double value(const double &density) override;

        Species getFilter() override { return _idSpecies;};

    private:
        // raw params
        Species _idSpecies;
        double _energyScale;
        double _lengthScale;
        double _density;

        // optimized for calculation
        double _totalPrefactor;
        double _inverseThetaSq;

        double derivative(const double &density);
    };

    REGISTER_PARSER("squared_exp", EamKernel, EamSE)
}

#endif
