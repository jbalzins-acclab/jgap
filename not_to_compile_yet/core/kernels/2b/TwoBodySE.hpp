#ifndef JGAP_TWOBODYSE_HPP
#define JGAP_TWOBODYSE_HPP

#include "TwoBodyKernel.hpp"
#include "../Kernel.hpp"
#include "core/cutoff/CutoffFunction.hpp"
#include "../../../data/atomic/AtomicStructure.hpp"
#include "../../../data/atomic/PredictionData.hpp"
#include "data/descriptors/kernels/TwoBodyIndex.hpp"
#include "io/parse/ParserRegistry.hpp"
#include "data/DataNode.hpp"

#include <string>
#include <optional>
#include <memory>

namespace jgap {

    class TwoBodySE : public TwoBodyKernel, Serializable {
    public:
        SETUP_PARSER_AND_SERIALIZATION(TwoBodyKernel, TwoBodySE, squared_exp_2b)
        ADD_PARSER(IKernel, TwoBodySE, squared_exp_2b)

        TwoBodySE(SpeciesPair species_pair, double energy_scale, double length_scale, double r, double f_cut,
                  std::optional<double> coeff = {});

        SpeciesPair getFilter() override { return id_pair_; }
        double crossCovariance(const std::shared_ptr<IKernel>& other) override;

    protected:
        double valueInternal(const double &r) const override;
        double derivativeInternal(const double &changing_r) const override;

    private:
        // raw params
        SpeciesPair id_pair_;
        double energy_scale_;
        double length_scale_;
        double r_;
        double descriptor_prefactors_; // essentially f_cut

        // optimized for calculation
        double total_prefactor_;
        double inverse_theta_sq_;
    };
}

#endif
