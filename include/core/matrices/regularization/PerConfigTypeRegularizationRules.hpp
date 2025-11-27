#ifndef PERCONFIGTYPEREGULARIZATIONRULES_HPP
#define PERCONFIGTYPEREGULARIZATIONRULES_HPP

#include <nlohmann/json.hpp>

#include "RegularizationRules.hpp"
#include "data/AtomicStructure.hpp"
#include "io/parse/ParserRegistry.hpp"

namespace jgap {
    class PerConfigTypeRegularizationRules : public RegularizationRules {
        public:
            static constexpr string TYPE = "per_config_type";

            PerConfigTypeRegularizationRules(const nlohmann::json &params);
            ~PerConfigTypeRegularizationRules() override = default;

            void fillSigmas(AtomicStructure &structure) override;

        private:
            double _energy;
            double _force;
            // Virial regularization split: isotropic (diagonals) vs anisotropic (off-diagonals)
            // If not provided, both default to the legacy single V value
            double _virialIso;
            double _virialAniso;

            // Order-sensitive list of rules. Supports both exact match (is) and substring (contains)
            // We prioritize exact matches over contains at evaluation time.
            struct Rule {
                bool exact;            // true => "is" rule, false => "contains" rule
                string key;            // value to compare with config_type
                double multiplier;     // multiplier to apply when matched
            };
            vector<Rule> _rules;
    };

    REGISTER_PARSER(RegularizationRules, PerConfigTypeRegularizationRules);
}

#endif