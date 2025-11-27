#ifndef JGAP_NULLPAIRFUNCTION_HPP
#define JGAP_NULLPAIRFUNCTION_HPP

#include "EamPairFunction.hpp"
#include "io/log/CurrentLogger.hpp"
#include "io/parse/ParserRegistry.hpp"

namespace jgap {
    class NullPairFunction : public EamPairFunction {
    public:
        static constexpr string TYPE = "null";

        NullPairFunction();
        NullPairFunction(const nlohmann::json params) {}

        double evaluate(double distance) override { return 0.0; }
        double differentiate(double distance) override { return 0.0; }
        string getType() override { return TYPE; }
        nlohmann::json serialize() override { return nlohmann::json(); }
    };
    REGISTER_PARSER(EamPairFunction, NullPairFunction)
}

#endif