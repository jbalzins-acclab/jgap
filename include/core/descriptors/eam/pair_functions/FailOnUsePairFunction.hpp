#ifndef JGAP_FAILONUSEPAIRFUNCTION_HPP
#define JGAP_FAILONUSEPAIRFUNCTION_HPP

#include "EamPairFunction.hpp"
#include "io/log/CurrentLogger.hpp"
#include "io/parse/ParserRegistry.hpp"

namespace jgap {
    class FailOnUsePairFunction : public EamPairFunction {
    public:
        static constexpr string TYPE = "fail_on_use";

        FailOnUsePairFunction();
        FailOnUsePairFunction(const nlohmann::json& json);

        double evaluate(double distance) override { doFail(); }
        double differentiate(double distance) override { doFail(); }
        string getType() override { return TYPE; }
        nlohmann::json serialize() override { return nlohmann::json(); }

    private:
        static void doFail() {
            JGAP_LOG_AND_THROW(
                "Attempted evaluation of a {} EAM pair function"
                "(try comparing species in data and what is specified in settings)",
                TYPE
                );
        }
    };
    REGISTER_PARSER(EamPairFunction, FailOnUsePairFunction)
}

#endif