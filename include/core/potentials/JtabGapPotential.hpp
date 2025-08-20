#ifndef JTABGAPPOTENTIAL_HPP
#define JTABGAPPOTENTIAL_HPP
#include "JgapPotential.hpp"

using namespace std;

namespace jgap {
    class JtabGapPotential : public Potential {
    public:
        JtabGapPotential(const nlohmann::json& params);
        ~JtabGapPotential() override = default;

        PotentialPrediction predict(const AtomicStructure &structure) override;

        nlohmann::json serialize() override;

        string getType() override {return "jtabgap";};
        double getCutoff() override;

        TabulationData tabulate(const TabulationParams &params) override;
    };
    REGISTER_PARSER("jtabgap", Potential, JtabGapPotential);
}

#endif //JTABGAPPOTENTIAL_HPP
