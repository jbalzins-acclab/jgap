#ifndef JGAP_FAILONUSEPAIRFUNCTION_HPP
#define JGAP_FAILONUSEPAIRFUNCTION_HPP

#include "EamPairFunction.hpp"
#include "io/log/CurrentLogger.hpp"
#include "io/parse/ParserRegistry.hpp"

namespace jgap {
    class FailOnUsePairFunction : public EamPairFunction, Serializable {
    public:
        SETUP_PARSER_AND_SERIALIZATION(EamPairFunction, FailOnUsePairFunction, fail_on_use)

        FailOnUsePairFunction();

        double evaluate(double distance) override { doFail(); }
        double differentiate(double distance) override { doFail(); }

    private:
        static void doFail() {
            JGAP_LOG_AND_THROW(
                "Attempted evaluation of a {} EAM std::pair function"
                "(try comparing species in data and what is specified in settings)",
                TYPE_ID
                );
        }
    };

    inline DataNode FailOnUsePairFunction::serialize() {
        return DataNode::object();
    }
}

#endif