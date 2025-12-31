#ifndef JGAP_PERCONFIGTYPEREGULARIZATIONRULES_HPP
#define JGAP_PERCONFIGTYPEREGULARIZATIONRULES_HPP

#include "RegularizationRules.hpp"
#include "../../../data/atomic/AtomicStructure.hpp"
#include "io/parse/ParserRegistry.hpp"

namespace jgap {
    class PerConfigTypeRegularizationRules : public RegularizationRules {
    public:
        SETUP_PARSER(RegularizationRules, PerConfigTypeRegularizationRules, multiplier_per_config_type);

        PerConfigTypeRegularizationRules(
            double energyPerAtom, double forceComponent, double virialIsoPerAtom, double virialAnisoPerAtom,
            const std::map<std::string, double> &exact, const std::map<std::string, double> &contains
            );
        ~PerConfigTypeRegularizationRules() override = default;

        void fillSigmas(std::vector<AtomicStructure>& structure) override;

    private:
        double _defaultEnergyPerAtom;
        double _defaultForceComponent;
        double _defaultVirialIsoPerAtom;
        double _defaultVirialAnisoPerAtom;

        std::map<std::string, double> _exact;
        std::map<std::string, double> _contains; // order-sensitive

        void fillSigmas(AtomicStructure &structure);
    };
}

#endif