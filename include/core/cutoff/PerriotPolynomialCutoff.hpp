#ifndef JGAP_PERRIOTPOLYNOMIALCUTOFF_HPP
#define JGAP_PERRIOTPOLYNOMIALCUTOFF_HPP

#include "CutoffFunction.hpp"
#include "io/parse/ParserRegistry.hpp"

namespace jgap {
    class PerriotPolynomialCutoff : public CutoffFunction {
    public:
        static constexpr std::string TYPE = "perriot";

        PerriotPolynomialCutoff(const DataNode& params);
        PerriotPolynomialCutoff(double rMin, double cutoff);

        std::string getType() override { return TYPE; }
        DataNode serialize() override;
        double getCutoff() override { return _cutoff; }

        double evaluate(double r) override;
        double differentiate(double r) override;

    private:
        double _cutoff;
        double _rMin;
        double _cutoffWidthInverse;
    };

    SETUP_PARSER(CutoffFunction, PerriotPolynomialCutoff)
}

#endif
