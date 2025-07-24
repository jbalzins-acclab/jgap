#ifndef PERRIOTPOLYNOMIALCUTOFF_HPP
#define PERRIOTPOLYNOMIALCUTOFF_HPP

#include "CutoffFunction.hpp"
#include "io/parse/ParserRegistry.hpp"

namespace jgap {
    class PerriotPolynomialCutoff : public CutoffFunction {
    public:
        explicit PerriotPolynomialCutoff(const nlohmann::json& params);
        explicit PerriotPolynomialCutoff(double rMin, double cutoff);

        string getType() override { return "perriot"; }
        nlohmann::json serialize() override;
        double getCutoff() override { return _cutoff; }

        double evaluate(double r) override;
        double differentiate(double r) override;

    private:
        double _cutoff;
        double _rMin;
        double _cutoffWidthInverse;
    };

    REGISTER_PARSER("perriot", CutoffFunction, PerriotPolynomialCutoff)
}

#endif //PERRIOTPOLYNOMIALCUTOFF_HPP
