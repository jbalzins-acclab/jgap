#ifndef ISOLATEDATOMPOTENTIAL_HPP
#define ISOLATEDATOMPOTENTIAL_HPP

#include <map>
#include <nlohmann/json.hpp>

#include "Potential.hpp"
#include "io/parse/ParserRegistry.hpp"


namespace jgap {
    class IsolatedAtomPotential : public Potential {
    public:
        static constexpr string TYPE = "isolated_atom";

        explicit IsolatedAtomPotential(const nlohmann::json& params);
        explicit IsolatedAtomPotential(const std::map<Species, double>& isolatedAtomEnergies, bool errorOnUnknown);

        nlohmann::json serialize() override;
        string getType() override { return TYPE; }
        CutoffRanges getCutoff() override { return {}; }

        Predictions predict(const AtomicStructure& structure) override;

        void tabulate(TabulationData& table) override;

    private:
        bool _errorOnUnknownSpecies;
        std::map<Species, double> _isolatedEnergies;
    };

    REGISTER_PARSER(Potential, IsolatedAtomPotential);
}

#endif //ISOLATEDATOMPOTENTIAL_HPP
