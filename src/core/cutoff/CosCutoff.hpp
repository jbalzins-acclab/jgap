#ifndef JGAP_DEFAULTCUTOFFFUNCTION_HPP
#define JGAP_DEFAULTCUTOFFFUNCTION_HPP

#include "CutoffFunction.hpp"
#include "io/parse/ParserRegistry.hpp"

namespace jgap {
    class CosCutoff : public CutoffFunction {
    public:
        CosCutoff(double cutoff, double cutoffTransitionWidth);
        explicit CosCutoff(const DataNode& params);

        static constexpr const char* TYPE = "coscutoff";

        std::string getType() override {return TYPE;};
        DataNode serialize() override;
        double getCutoff() override { return _cutoff; }

        ~CosCutoff() override = default;

        double evaluate(double r) override;
        double differentiate(double r) override;

    private:
        double _cutoff;
        double _cutoffTransitionWidth;
        double _cutoffTransitionWidthInverse;
    };

    SETUP_PARSER(CutoffFunction, CosCutoff)
}
#endif
