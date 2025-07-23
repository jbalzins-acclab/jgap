//
// Created by Jegors Balzins on 21.7.2025.
//

#ifndef FORCESCALESIGMARULES_HPP
#define FORCESCALESIGMARULES_HPP
#include "SigmaRules.hpp"
#include "io/parse/ParserRegistry.hpp"

namespace jgap {
    class ForceScaleSigmaRules : public SigmaRules {

    public:
        ForceScaleSigmaRules(const nlohmann::json &params);
        ~ForceScaleSigmaRules() override = default;

        void fillSigmas(AtomicStructure &structure) override;

    private:
        double _EFratio;
        double _VFratio;
        double _sigmaFmin;
        double _scalingFactor;
        double _referenceF;
        double _referenceSigmaF;
        bool _componentWise;

        double calculateSigma(double F) const;
    };

    REGISTER_PARSER("force_scale", SigmaRules, ForceScaleSigmaRules);
}



#endif //SCALEDFORCESIGMARULES_HPP
