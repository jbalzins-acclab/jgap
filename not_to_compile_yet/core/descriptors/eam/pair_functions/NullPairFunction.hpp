#ifndef JGAP_NULLPAIRFUNCTION_HPP
#define JGAP_NULLPAIRFUNCTION_HPP

#include "EamPairFunction.hpp"
#include "io/log/CurrentLogger.hpp"
#include "io/parse/ParserRegistry.hpp"

namespace jgap {
    class NullPairFunction : public EamPairFunction, Serializable {
    public:
        SETUP_PARSER_AND_SERIALIZATION(EamPairFunction, NullPairFunction, null)

        NullPairFunction() = default;

        double evaluate(double distance) override { return 0.0; }
        double differentiate(double distance) override { return 0.0; }
    };

    inline std::shared_ptr<NullPairFunction> NullPairFunction::fromDataNode(const DataNode &params) {
        return std::make_shared<NullPairFunction>();
    }

    inline DataNode NullPairFunction::serialize() { return {}; }
}

#endif