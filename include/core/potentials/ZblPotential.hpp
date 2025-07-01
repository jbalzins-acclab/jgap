//
// Created by Jegors Balzins on 22.6.2025.
//

#ifndef ZBLPOTENTIAL_HPP
#define ZBLPOTENTIAL_HPP

#include <nlohmann/json.hpp>

#include "Potential.hpp"
#include "core/cutoff/CutoffFunction.hpp"
#include "io/log/StdoutLogger.hpp"
#include "utils/AtomicNumbers.hpp"
#include "utils/Utils.hpp"

namespace jgap {
    class ZblPotential : public Potential {
    public:
        explicit ZblPotential(nlohmann::json dmolFitCoefficients, double cutoffLow = 1.2, double cutoffHigh = 2.2);
        ~ZblPotential() override = default;

        PotentialPrediction predict(const AtomicStructure& structure) override;
        nlohmann::json serialize() override;

    private:
        map<SpeciesPair, array<double, 6>> _dmolFitCoefficients;
        shared_ptr<CutoffFunction> _cutoffFunction;
        double _cutoff;

        const double _eps = 8.854187817e-12;
        const double _electronCharge = 1.60217657e-19;

        double zbl_eV(const SpeciesPair& speciesPair, double r);
        double zblWithCutoff_eV(const SpeciesPair& speciesPair, double r);
        double zblWithCutoffDerivative_eV_per_Ang(const SpeciesPair& speciesPair, double r);
    };
}

#endif //ZBLPOTENTIAL_HPP
