#ifndef JGAP_DEFAULTCUTOFFFUNCTION_HPP
#define JGAP_DEFAULTCUTOFFFUNCTION_HPP

#include "CutoffFunction.hpp"
#include "io/parse/ParserRegistry.hpp"

namespace jgap {
    class CosCutoff : public CutoffFunction {
    public:
        static constexpr string TYPE = "coscutoff";

        CosCutoff(nlohmann::json params);
        CosCutoff(double cutoff, double rMin);

        string getType() override {return TYPE;};
        nlohmann::json serialize() override;
        double getCutoff() override { return _cutoff; }

        ~CosCutoff() override = default;

        double evaluate(double r) override;
        double differentiate(double r) override;

    private:
        double _cutoff;
        double _cutoffTransitionWidth;
        double _cutoffTransitionWidthInverse;
    };

    REGISTER_PARSER(CutoffFunction, CosCutoff)
}
#endif
