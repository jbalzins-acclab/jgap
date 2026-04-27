#ifndef THREEBODYSE_HPP
#define THREEBODYSE_HPP

#include "core/cutoff/CutoffFunction.hpp"
#include "core/descriptors/Kernel.hpp"
#include "data/descriptors/kernels/ThreeBodyIndex.hpp"
#include "io/parse/ParserRegistry.hpp"

#include <nlohmann/json.hpp>

#include "io/parse/ParserRegistry.hpp"

namespace jgap {

    using ThreeBodyKernel = Kernel<SpeciesTriplet, ThreeBodyIndex, ThreeBodyDescriptorData>;

    class ThreeBodySE : public ThreeBodyKernel {
    public:
        ThreeBodySE(SpeciesTriplet idTriplet, double energyScale, Vector3 lengthScales, Vector3 q, double fCut,
                    optional<double> coeff = {});

        ThreeBodySE(const nlohmann::json& params);
        string getType() override { return "squared_exp"; }
        nlohmann::json serialize() override;

        double crossCovariance(const shared_ptr<IKernel> &other) override;
        Covariance covariance(const AtomicStructure &structure, const ThreeBodyIndex &index) override;
        double value(const ThreeBodyDescriptorData &t) override;

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

        [[nodiscard]] double valueInternal(const Vector3 &q) const;
        [[nodiscard]] Vector3 gradientInternal(const Vector3 &q) const;
    };

    REGISTER_PARSER("squared_exp", ThreeBodyKernel, ThreeBodySE);
}

#endif //THREEBODYSE_HPP
