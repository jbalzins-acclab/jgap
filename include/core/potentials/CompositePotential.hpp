#ifndef COMPOSITEPOTENTIAL_HPP
#define COMPOSITEPOTENTIAL_HPP

#include "Potential.hpp"
#include "io/parse/ParserRegistry.hpp"

namespace jgap {
    class CompositePotential : public Potential {
    public:
        ~CompositePotential() override = default;

        explicit CompositePotential(map<string, shared_ptr<Potential>> &potentials);
        explicit CompositePotential(const nlohmann::json& params);
        nlohmann::json serialize() override;
        string getType() override { return "composite"; };
        double getCutoff() override;

        PotentialPrediction predict(const AtomicStructure &structure) override;

        TabulationData tabulate(const TabulationParams& params) override;

    private:
        map<string, shared_ptr<Potential>> _potentials;
    };

    REGISTER_PARSER("composite", Potential, CompositePotential)
}

#endif //COMPOSITEPOTENTIAL_HPP
