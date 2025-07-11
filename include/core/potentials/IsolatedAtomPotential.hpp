//
// Created by Jegors Balzins on 22.6.2025.
//

#ifndef ISOLATEDATOMPOTENTIAL_HPP
#define ISOLATEDATOMPOTENTIAL_HPP

#include <map>
#include <nlohmann/json.hpp>

#include "Potential.hpp"
#include "io/parse/ParserRegistry.hpp"

using namespace std;

namespace jgap {
    class IsolatedAtomPotential : public Potential {
    public:
        ~IsolatedAtomPotential() override = default;

        explicit IsolatedAtomPotential(const nlohmann::json& params);
        nlohmann::json serialize() override;
        string getType() override { return "isolated_atom"; }
        double getCutoff() override { return 0.0; }

        PotentialPrediction predict(const AtomicStructure& structure) override;


    private:
        bool _errorOnUnknownSpecies;
        map<Species, double> _isolatedEnergies;
    };

    REGISTER_PARSER("isolated_atom", Potential, IsolatedAtomPotential);
}

#endif //ISOLATEDATOMPOTENTIAL_HPP
