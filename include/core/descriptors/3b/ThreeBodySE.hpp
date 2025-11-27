#ifndef JGAP_THREEBODYSE_HPP
#define JGAP_THREEBODYSE_HPP

#include "core/cutoff/CutoffFunction.hpp"
#include "core/descriptors/Kernel.hpp"
#include "data/descriptors/kernels/ThreeBodyIndex.hpp"
#include "io/parse/ParserRegistry.hpp"

#include <nlohmann/json.hpp>

#include "ThreeBodyKernel.hpp"
#include "io/parse/ParserRegistry.hpp"

namespace jgap {
    class ThreeBodySE : public ThreeBodyKernel {
    public:
        static constexpr std::string TYPE = "squared_exp";
        ThreeBodySE(SpeciesTriplet idTriplet, double energyScale, Vector3 lengthScales, Vector3 q, double fCut);

        ThreeBodySE(const nlohmann::json& params);
        std::string getType() override { return TYPE; }
        nlohmann::json serialize() override;

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

    REGISTER_PARSER(ThreeBodyKernel, ThreeBodySE);
}

#endif
