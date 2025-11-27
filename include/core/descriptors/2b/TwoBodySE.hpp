#ifndef JGAP_TWOBODYSE_HPP
#define JGAP_TWOBODYSE_HPP

#include "TwoBodyKernel.hpp"
#include "../Kernel.hpp"
#include "core/cutoff/CutoffFunction.hpp"
#include "data/AtomicStructure.hpp"
#include "data/PredictionData.hpp"
#include "data/descriptors/kernels/TwoBodyIndex.hpp"
#include "io/parse/ParserRegistry.hpp"

namespace jgap {


    class TwoBodySE : public TwoBodyKernel {
    public:
        static constexpr string TYPE = "squared_exp";
        TwoBodySE(SpeciesPair speciesPair, double energyScale, double lengthScale, double r, double fCut,
                  optional<double> coeff = {});
        TwoBodySE(const nlohmann::json &params);
        string getType() override { return TYPE; }
        nlohmann::json serialize() override;

        SpeciesPair getFilter() override { return _idPair; }

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

    REGISTER_PARSER(TwoBodyKernel, TwoBodySE)
}

#endif
