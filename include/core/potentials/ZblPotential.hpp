//
// Created by Jegors Balzins on 22.6.2025.
//

#ifndef ZBLPOTENTIAL_HPP
#define ZBLPOTENTIAL_HPP

#include <nlohmann/json.hpp>

#include "Potential.hpp"
#include "core/cutoff/CutoffFunction.hpp"
#include "io/parse/ParserRegistry.hpp"
#include "utils/AtomicNumbers.hpp"
#include "utils/Utils.hpp"

#define DEFAULT_ZBL_CUTOFF 2.2
#define DEFAULT_ZBL_RMIN 1.2

namespace jgap {
    class ZblPotential : public Potential {
    public:
        explicit ZblPotential(const nlohmann::json& zblParams);
        ~ZblPotential() override = default;

        PotentialPrediction predict(const AtomicStructure& structure) override;
        nlohmann::json serialize() override;
        string getType() override { return "zbl"; }
        double getCutoff() override { return _cutoff; }

        TabulationData tabulate(const TabulationParams &params) override;

    private:
        double _cutoff;
        string _dmolFile;
        map<SpeciesPair, array<double, 6>> _dmolFitCoefficients;
        shared_ptr<CutoffFunction> _cutoffFunction;

        const double _eps = 8.854187817e-12;
        const double _electronCharge = 1.60217657e-19;

        static string getDefaultDmolFilePath();

        double zbl_eV(const SpeciesPair& speciesPair, double r);
        double zblWithCutoff_eV(const SpeciesPair& speciesPair, double r);
        double zblWithCutoffDerivative_eV_per_Ang(const SpeciesPair& speciesPair, double r);

        void parseDmolFitCoefficients();
    };

    REGISTER_PARSER("zbl", Potential, ZblPotential)
}

#endif //ZBLPOTENTIAL_HPP
