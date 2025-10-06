#ifndef TWOBODYSE_HPP
#define TWOBODYSE_HPP

#include "../Kernel.hpp"
#include "core/cutoff/CutoffFunction.hpp"
#include "data/descriptors/kernels/TwoBodyIndex.hpp"
#include "io/parse/ParserRegistry.hpp"

namespace jgap {

    using TwoBodyKernel = Kernel<SpeciesPair, TwoBodyIndex, TwoBodyDescriptorData>;

    class TwoBodySE : public TwoBodyKernel {
    public:
        TwoBodySE(SpeciesPair speciesPair, double energyScale, double lengthScale, double r, double fCut,
                  optional<double> coeff = {});
        TwoBodySE(const nlohmann::json &params);
        string getType() override { return "squared_exp"; }
        nlohmann::json serialize() override;

        Covariance covariance(const AtomicStructure &structure, const TwoBodyIndex &index) override;
        double value(const TwoBodyDescriptorData &r) override;
        SpeciesPair getFilter() override { return _idPair;};

        double crossCovariance(const shared_ptr<IKernel>& other) override;

    private:
        // raw params
        SpeciesPair _idPair;
        double _energyScale;
        double _lengthScale;
        double _r;
        double _descriptorPrefactors; // essentially f_cut

        // optimized for calculation
        double _totalPrefactor;
        double _inverseThetaSq;

        double valueInternal(const double &r) const;
        double derivativeInternal(const double &changingR) const;
    };


    REGISTER_PARSER("squared_exp", TwoBodyKernel, TwoBodySE)
}

#endif
