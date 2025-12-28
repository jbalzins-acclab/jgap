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
        SETUP_PARSER_AND_SERIALIZATION(TwoBodyKernel, TwoBodySE, squared_exp)

        TwoBodySE(SpeciesPair speciesPair, double energyScale, double lengthScale, double r, double fCut,
                  std::optional<double> coeff = {});

        SpeciesPair getFilter() override { return _idPair; }
        double crossCovariance(const std::shared_ptr<IKernel>& other) override;

    protected:
        double valueInternal(const double &r) const override;
        double derivativeInternal(const double &changingR) const override;

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
    };
}

#endif
