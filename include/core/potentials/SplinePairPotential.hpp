#ifndef SPLINEPAIRPOTENTIAL_HPP
#define SPLINEPAIRPOTENTIAL_HPP

#include "Potential.hpp"
#include "io/parse/ParserRegistry.hpp"

namespace jgap {
    class SplinePairPotential : public Potential {
    public:
        class NaturalCubicSpline {
        // NOTE: ChatGPT implementation - TODO: check
        public:
            NaturalCubicSpline(nlohmann::json params);
            NaturalCubicSpline(vector<double> r, vector<double> E);

            double evaluate(double r) const;
            double derivative(double r) const;
            double getCutoff() const { return _r.back(); };

            nlohmann::json serialize() const;

        private:
            vector<double> _r, _a, _b, _c, _d;

            void init(const vector<double> &r, const vector<double> &E);
            size_t findInterval(double r) const;
        };

        explicit SplinePairPotential(nlohmann::json params);
        explicit SplinePairPotential(map<SpeciesPair, pair<vector<double>, vector<double>> > points);
        ~SplinePairPotential() override = default;

        PotentialPrediction predict(const AtomicStructure &structure) override;

        nlohmann::json serialize() override;
        string getType() override { return "spline_pairpot"; }
        double getCutoff() override;

        TabulationData tabulate(const TabulationParams &params) override;

    private:
        map<SpeciesPair, shared_ptr<NaturalCubicSpline>> _perSpeciesInterpolators;
    };
    REGISTER_PARSER("spline_pairpot", Potential, SplinePairPotential)
}

#endif
