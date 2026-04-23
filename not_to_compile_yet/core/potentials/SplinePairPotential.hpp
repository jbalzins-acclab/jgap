#ifndef SPLINEPAIRPOTENTIAL_HPP
#define SPLINEPAIRPOTENTIAL_HPP

#include "Potential.hpp"
#include "../../data/atomic/AtomicStructure.hpp"
#include "core/tabulation/Tabulatable.hpp"
#include "io/Serializable.hpp"
#include "io/parse/ParserRegistry.hpp"

namespace jgap {
    class SplinePairPotential : public Potential, Serializable, Tabulatable {
    public:
        SETUP_PARSER_AND_SERIALIZATION(Potential, SplinePairPotential, spline_pairpot)

        class NaturalCubicSpline {
        public:
            NaturalCubicSpline(const DataNode& params);
            NaturalCubicSpline(const std::vector<double>& r, const std::vector<double>& E);

            double evaluate(double r) const;
            double derivative(double r) const;
            double getCutoff() const { return r_.back(); };

            DataNode serialize() const;

        private:
            std::vector<double> r_, energies_, b_, c_, d_;

            void init(const std::vector<double> &r, const std::vector<double> &E);
            std::size_t findInterval(double r) const;
        };

        SplinePairPotential(const std::map<SpeciesPair, std::pair<std::vector<double>, std::vector<double>> >& points);
        ~SplinePairPotential() override = default;

        Predictions predict(const AtomicStructure &structure) override;

        CutoffRanges getCutoff() override;

        void tabulate(TabulationData& table) override;

    private:
        std::map<SpeciesPair, std::shared_ptr<NaturalCubicSpline>> per_species_interpolators_;
    };
}

#endif
