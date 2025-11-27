#ifndef ZBLPOTENTIAL_HPP
#define ZBLPOTENTIAL_HPP

#include <nlohmann/json.hpp>
#include <set>

#include "Potential.hpp"
#include "core/cutoff/CutoffFunction.hpp"
#include "io/parse/ParserRegistry.hpp"
#include "utils/AtomicNumbers.hpp"
#include "utils/Utils.hpp"

#define DEFAULT_ZBL_CUTOFF 2.2
#define DEFAULT_ZBL_R_MIN 1.2

namespace jgap {
    class ZblPotential : public Potential {
    public:
        static constexpr string TYPE = "zbl";

        ZblPotential(const nlohmann::json& zblParams);
        ~ZblPotential() override = default;

        Predictions predict(const AtomicStructure& structure) override;
        nlohmann::json serialize() override;
        string getType() override { return TYPE; }
        CutoffRanges getCutoff() override { return CutoffRanges{.twoBody = _cutoff}; }

        void tabulate(TabulationData& table) override;

    private:
        static constexpr double _MINIMAL_TABULATED_R = 1e-4;
        static constexpr double _EPSILON = 8.854187817e-12;
        static constexpr double _ELECTRON_CHARGE = 1.60217657e-19;

        double _cutoff;
        string _coeffFileName;

        set<Species> _relevantSpecies{};

        std::map<SpeciesPair, array<double, 6>> _dmolFitCoefficients;
        shared_ptr<CutoffFunction> _cutoffFunction;

        static string getResourcesCoeffFilePath(const string& fileName = "dmol-fit.json");

        double zbl_eV(const SpeciesPair& speciesPair, double r);
        double zblWithCutoff_eV(const SpeciesPair& speciesPair, double r);
        double zblWithCutoffDerivative_eV_per_Ang(const SpeciesPair& speciesPair, double r);
    };

    REGISTER_PARSER(Potential, ZblPotential)
}

#endif //ZBLPOTENTIAL_HPP
