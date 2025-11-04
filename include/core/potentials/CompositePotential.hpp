#ifndef JGAP_COMPOSITEPOTENTIAL_HPP
#define JGAP_COMPOSITEPOTENTIAL_HPP

#include "Potential.hpp"
#include "data/AtomicStructure.hpp"
#include "data/PredictionData.hpp"
#include "io/parse/ParserRegistry.hpp"

namespace jgap {

    class CompositePotential : public Potential {
    public:
        static constexpr string TYPE = "composite";

        ~CompositePotential() override = default;

        explicit CompositePotential(const map<string, shared_ptr<Potential>>& potentials);
        explicit CompositePotential(const nlohmann::json& params);
        nlohmann::json serialize() override;
        string getType() override { return TYPE; }
        CutoffRanges getCutoff() override;

        Predictions predict(const AtomicStructure &structure) override;

        void tabulate(TabulationData& table) override;

    private:
        map<string, shared_ptr<Potential>> _potentials;
    };

    REGISTER_PARSER(Potential, CompositePotential)
}

#endif //COMPOSITEPOTENTIAL_HPP
