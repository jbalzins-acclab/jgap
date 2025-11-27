#ifndef SIMPLEREGULARIZATIONRULES_HPP
#define SIMPLEREGULARIZATIONRULES_HPP

#include "RegularizationRules.hpp"
#include "io/parse/ParserRegistry.hpp"

namespace jgap {
    class SimpleRegularizationRules : public RegularizationRules {
    public:
        static constexpr string TYPE = "simple";

        SimpleRegularizationRules(const nlohmann::json &params);
        ~SimpleRegularizationRules() override = default;

        void fillSigmas(AtomicStructure &structure) override;

    private:
        double _defaultEPerAtom;
        double _defaultF;
        double _defaultVirialsPerAtom;
        double _liquidMultiplier;
        double _shortRangeMultiplier;
    };

    REGISTER_PARSER(RegularizationRules, SimpleRegularizationRules);
}
#endif
