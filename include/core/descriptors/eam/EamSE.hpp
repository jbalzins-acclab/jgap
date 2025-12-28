#ifndef EAMKERNEL_HPP
#define EAMKERNEL_HPP

#include "EamKernel.hpp"
#include "../Kernel.hpp"
#include "data/descriptors/kernels/EamKernelIndex.hpp"
#include "io/parse/ParserRegistry.hpp"

namespace jgap {
    class EamSE : public EamKernel, Serializable {
    public:
        SETUP_PARSER_AND_SERIALIZATION(EamKernel, EamSE, squared_exp)

        EamSE(Species species, double energyScale, double lengthScale, double density,
              std::optional<double> coeff = {});

        double crossCovariance(const std::shared_ptr<IKernel> &other) override;

        Species getFilter() override { return id_species_; }
        std::pair<double, double> getDensityRange() override { return {density_, density_}; }

    private:
        // raw params
        Species id_species_;
        double energy_scale_;
        double length_scale_;
        double density_;

        // optimized for calculation
        double total_prefactor_;
        double inverse_theta_sq;

        double valueInternal(const double &density) const override;
        double derivativeInternal(const double &density) const override;
    };
}

#endif
