//
// Created by Jegors Balzins on 22.6.2025.
//

#ifndef ISOLATEDATOMPOTENTIAL_HPP
#define ISOLATEDATOMPOTENTIAL_HPP

#include <map>
#include <nlohmann/json.hpp>

#include "Potential.hpp"

using namespace std;

namespace jgap {
    class IsolatedAtomPotential : public Potential {
    public:
        explicit IsolatedAtomPotential(const vector<AtomicStructure> &containsIsolatedAtoms,
                                       bool errorOnUnknownSpecies = true);
        ~IsolatedAtomPotential() override = default;

        PotentialPrediction predict(const AtomicStructure& structure) override;

        nlohmann::json serialize() override;

    private:
        bool _errorOnUnknownSpecies;
        map<Species, double> _isolatedEnergies;
    };
}

#endif //ISOLATEDATOMPOTENTIAL_HPP
