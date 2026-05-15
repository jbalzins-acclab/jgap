#ifndef JGAP_PERRIOTPOLYNOMIALCUTOFF_HPP
#define JGAP_PERRIOTPOLYNOMIALCUTOFF_HPP

#include "CutoffFunction.hpp"
#include "core/Real.hpp"
#include "io/parse/ParserRegistry.hpp"

namespace jgap {
    class PerriotPolynomialCutoff : public CutoffFunction {
    public:
        static constexpr std::string TYPE = "perriot";

        PerriotPolynomialCutoff(Real r_min, Real cutoff);

        Real getCutoff() const override { return cutoff; }

        Real evaluate(Real r) const override;
        Real differentiate(Real r) const override;

    private:
        Real cutoff;
        Real r_min;
        Real cutoff_width_inverse;
    };

    //SETUP_PARSER(CutoffFunction, PerriotPolynomialCutoff)
}

#endif
