#ifndef JGAP_PERCONFIGTYPEREGULARIZATIONRULES_HPP
#define JGAP_PERCONFIGTYPEREGULARIZATIONRULES_HPP

#include <nlohmann/json.hpp>

#include "RegularizationRules.hpp"
#include "io/parse/ParserRegistry.hpp"

namespace jgap {
    class PerConfigTypeRegularizationRules : public RegularizationRules {
        public:
            explicit PerConfigTypeRegularizationRules(const nlohmann::json &params);
            ~PerConfigTypeRegularizationRules() override = default;

            void fillSigmas(AtomicStructure &structure) override;

        private:
            double _E;
            double _F;
            double _V;
            vector<pair<string, double>> _multipliersPerKeyWord;
    };

    REGISTER_PARSER("per_ct", RegularizationRules, PerConfigTypeRegularizationRules);
}

#endif