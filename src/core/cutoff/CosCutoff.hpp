#ifndef JGAP_DEFAULTCUTOFFFUNCTION_HPP
#define JGAP_DEFAULTCUTOFFFUNCTION_HPP

#include "CutoffFunction.hpp"
#include "core/Real.hpp"
#include "io/parse/ParserRegistry.hpp"

namespace jgap {
    class CosCutoff : public CutoffFunction {
    public:
        CosCutoff(Real cutoff, Real cutoff_transition_width);
        static constexpr const char* TYPE = "coscutoff";

        Real getCutoff() const override { return cutoff; }

        ~CosCutoff() override = default;

        Real evaluate(Real r) const override;
        Real differentiate(Real r) const override;

    private:
        Real cutoff;
        Real cutoff_transition_width;

        Real cutoff_transition_width_inverse;
    };

    //SETUP_PARSER(CutoffFunction, CosCutoff)
}
#endif
