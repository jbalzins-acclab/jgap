#ifndef JGAP_THREEBODYSE_HPP
#define JGAP_THREEBODYSE_HPP

#include "core/cutoff/CutoffFunction.hpp"
#include "../Kernel.hpp"
#include "data/descriptors/kernels/ThreeBodyIndex.hpp"
#include "io/parse/ParserRegistry.hpp"

#include "ThreeBodyKernel.hpp"
#include "io/parse/ParserRegistry.hpp"

namespace jgap {
    class ThreeBodySE : public ThreeBodyKernel, Serializable {
    public:
        SETUP_PARSER_AND_SERIALIZATION(ThreeBodyKernel, ThreeBodySE, squared_exp);

        ThreeBodySE(SpeciesTriplet idTriplet, double energyScale, Vector3 lengthScales, Vector3 q, double fCut);

        double crossCovariance(const std::shared_ptr<IKernel> &other) override;

        SpeciesTriplet getFilter() override { return _idTriplet; }

    private:
        // raw params
        SpeciesTriplet _idTriplet;
        double _energyScale;
        Vector3 _lengthScale;
        Vector3 _q;
        double _descriptorPrefactors; // essentially f_cut

        // optimized for calculation
        double _totalPrefactor;
        Vector3 _inverseThetaSq;

        double valueInternal(const Vector3 &q) const;
        Vector3 gradientInternal(const Vector3 &q) const;
    };

    ;
}

#endif
