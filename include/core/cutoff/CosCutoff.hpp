//
// Created by Jegors Balzins on 9.7.2025.
//

#ifndef DEFAULTCUTOFFFUNCTION_HPP
#define DEFAULTCUTOFFFUNCTION_HPP

#include "CutoffFunction.hpp"
#include "io/parse/ParserRegistry.hpp"

namespace jgap {
    class CosCutoff : public CutoffFunction {
    public:

        explicit CosCutoff(nlohmann::json params);

        string getType() override {return "coscutoff";};
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

    REGISTER_PARSER("coscutoff", CutoffFunction, CosCutoff)
}
#endif //DEFAULTCUTOFFFUNCTION_HPP
