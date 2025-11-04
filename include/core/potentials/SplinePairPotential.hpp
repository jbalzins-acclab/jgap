#ifndef SPLINEPAIRPOTENTIAL_HPP
#define SPLINEPAIRPOTENTIAL_HPP

#include "Potential.hpp"
#include "data/AtomicStructure.hpp"
#include "data/PredictionData.hpp"
#include "io/parse/ParserRegistry.hpp"

namespace jgap {
    class SplinePairPotential : public Potential {
    public:
        static constexpr string TYPE = "spline_pairpot";

        class NaturalCubicSpline {
        public:
            NaturalCubicSpline(nlohmann::json params);
            NaturalCubicSpline(const vector<double>& r, const vector<double>& E);

            double evaluate(double r) const;
            double derivative(double r) const;
            double getCutoff() const { return _r.back(); };

            nlohmann::json serialize() const;

        private:
            vector<double> _r, _energies, _b, _c, _d;

            void init(const vector<double> &r, const vector<double> &E);
            size_t findInterval(double r) const;
        };

        SplinePairPotential(const nlohmann::json& params);
        SplinePairPotential(const map<SpeciesPair, pair<vector<double>, vector<double>> >& points);
        ~SplinePairPotential() override = default;

        Predictions predict(const AtomicStructure &structure) override;

        nlohmann::json serialize() override;
        string getType() override { return TYPE; }
        CutoffRanges getCutoff() override;

        void tabulate(TabulationData& table) override;

    private:
        map<SpeciesPair, shared_ptr<NaturalCubicSpline>> _perSpeciesInterpolators;
    };

    REGISTER_PARSER(Potential, SplinePairPotential)
}

#endif
