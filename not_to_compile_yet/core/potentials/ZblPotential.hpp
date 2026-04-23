#ifndef ZBLPOTENTIAL_HPP
#define ZBLPOTENTIAL_HPP

#include <set>

#include "Potential.hpp"
#include "core/cutoff/CutoffFunction.hpp"
#include "core/tabulation/Tabulatable.hpp"
#include "io/Serializable.hpp"
#include "io/parse/ParserRegistry.hpp"
#include "utils/AtomicNumbers.hpp"
#include "utils/Utils.hpp"

#define DEFAULT_ZBL_CUTOFF 2.2
#define DEFAULT_ZBL_R_MIN 1.2

namespace jgap {
    class ZblPotential : public Potential, Serializable, Tabulatable {
    public:
        SETUP_PARSER_AND_SERIALIZATION(Potential, ZblPotential, zbl)

        ~ZblPotential() override = default;
        ZblPotential(double cutoff, std::string coeff_file_name, std::set<Species> relevant_species);

        Predictions predict(const AtomicStructure& structure) override;
        CutoffRanges getCutoff() override { return CutoffRanges{.twoBody = cutoff_}; }

        void tabulate(TabulationData& table) override;

    private:
        static constexpr double MINIMAL_TABULATED_R_ = 1e-4;
        static constexpr double EPSILON_ = 8.854187817e-12;
        static constexpr double ELECTRON_CHARGE_ = 1.60217657e-19;

        double cutoff_;
        std::string coeff_file_name_;

        std::set<Species> relevant_species_{};

        std::map<SpeciesPair, std::array<double, 6>> dmol_fit_coefficients_;
        std::shared_ptr<CutoffFunction> cutoff_function_;

        static std::string getResourcesCoeffFilePath(const std::string& fileName = "dmol-fit.json");

        double zbl_eV(const SpeciesPair& speciesPair, double r);
        double zblWithCutoff_eV(const SpeciesPair& speciesPair, double r);
        double zblWithCutoffDerivative_eV_per_Ang(const SpeciesPair& speciesPair, double r);
    };
}

#endif
