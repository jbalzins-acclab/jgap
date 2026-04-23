#ifndef JGAP_SIMPLEREGULARIZATIONRULES_HPP
#define JGAP_SIMPLEREGULARIZATIONRULES_HPP

#include "RegularizationRules.hpp"
#include "io/parse/ParserRegistry.hpp"

namespace jgap {
    class SimpleRegularizationRules : public RegularizationRules {
    public:
        SETUP_PARSER(RegularizationRules, SimpleRegularizationRules, simple);

        ~SimpleRegularizationRules() override = default;
        SimpleRegularizationRules(double energyPerAtom, double forceComponent,
                                  double virialIsoPerAtom, double virialAnisoPerAtom,
                                  double liquidMultiplier, double shortRangeMultiplier);

        void fillSigmas(std::vector<AtomicStructure> &structures) override;

    private:
        double _defaultEnergyPerAtom;
        double _defaultForceComponent;
        double _defaultVirialIsoPerAtom;
        double _defaultVirialAnisoPerAtom;

        double _liquidMultiplier;
        double _shortRangeMultiplier;

        void fillSigmas(AtomicStructure &structure);
    };
}
#endif
