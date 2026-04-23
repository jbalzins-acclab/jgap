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

        double getCutoff() override { return cutoff_; }

        double evaluate(double r) override;
        double differentiate(double r) override;

    private:
        double cutoff_;
        double r_min_;
        double cutoff_width_inverse_;
    };

    //SETUP_PARSER(CutoffFunction, PerriotPolynomialCutoff)
}

#endif
