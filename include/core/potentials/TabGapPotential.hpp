#ifndef JGAP_TABGAPPOTENTIAL_HPP
#define JGAP_TABGAPPOTENTIAL_HPP

#include <io/parse/ParserRegistry.hpp>

#include "Potential.hpp"
#include "data/AtomicStructure.hpp"
#include "data/PredictionData.hpp"

namespace jgap {

    class TabGapPotential : public Potential {
    public:
        static constexpr string TYPE = "tabgap";

        ~TabGapPotential() override = default;
        TabGapPotential(const nlohmann::json& params);
        TabGapPotential(TabulationData tabulationData, optional<string> );

        Predictions predict(const AtomicStructure &structure) override;

        nlohmann::json serialize() override;

        string getType() override { return TYPE; }
        CutoffRanges getCutoff() override;

        void tabulate(TabulationData& table) override;

        void saveFiles(string fileNamePrefix);

    private:
        friend class TabGapIO;

        map<Species, double> _isolatedEnergies;

        map<SpeciesPair, Grid1d> _splineCoeffs2b;
        map<SpeciesTriplet, Grid3d> _splineCoeffs3b;

        map<Species, Grid1d> _splineCoeffsEam;
        map<ContributorReceiverSpecies, Grid1d> _splineCoeffsEamPf;

        double eval2b(double r);
        double eval3b(double r_ij, double r_ik, double cos_theta);
        double evalEam(double density);
    };

    REGISTER_PARSER(Potential, TabGapPotential)
}

#endif
