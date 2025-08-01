#ifndef SIMPLEREGULARIZATIONRULES_HPP
#define SIMPLEREGULARIZATIONRULES_HPP

#include "RegularizationRules.hpp"
#include "io/parse/ParserRegistry.hpp"

namespace jgap {
    class SimpleRegularizationRules : public RegularizationRules {
    public:
        SimpleRegularizationRules(const nlohmann::json &params);

        ~SimpleRegularizationRules() override = default;

        void fillSigmas(AtomicStructure &structure) override;

    private:
        double _defaultEPerAtom;
        double _defaultF;
        double _defaultVirials;
        double _liquidMultiplier;
        double _shortRangeMultiplier;
    };

    REGISTER_PARSER("simple", RegularizationRules, SimpleRegularizationRules);
}
#endif
