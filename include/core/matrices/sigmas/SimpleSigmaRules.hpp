#ifndef SIMPLESIGMARULES_HPP
#define SIMPLESIGMARULES_HPP

#include "SigmaRules.hpp"
#include "io/parse/ParserRegistry.hpp"

namespace jgap {
    class SimpleSigmaRules : public SigmaRules {
    public:
        SimpleSigmaRules(const nlohmann::json &params);

        ~SimpleSigmaRules() override = default;

        void fillSigmas(AtomicStructure &structure) override;

    private:
        double _defaultEPerAtom;
        double _defaultF;
        double _defaultVirials;
        double _liquidMultiplier;
        double _shortRangeMultiplier;
    };

    REGISTER_PARSER("simple", SigmaRules, SimpleSigmaRules);
}
#endif //SIMPLESIGMARULES_HPP
