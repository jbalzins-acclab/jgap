#ifndef FSGENPAIRFUNCTION_HPP
#define FSGENPAIRFUNCTION_HPP

#include "EamPairFunction.hpp"
#include "io/parse/ParserRegistry.hpp"
#include "data/DataNode.hpp"

#include <string>

namespace jgap {
    class FSGenPairFunction : public EamPairFunction, Serializable {
    public:
        SETUP_PARSER_AND_SERIALIZATION(EamPairFunction, FSGenPairFunction, fsgen)

        FSGenPairFunction(const double cutoff, const double degree, const double prefactor = 1.0)
            : EamPairFunction(cutoff, prefactor), degree_(degree) {
            cutoff_inverse_ = 1.0 / cutoff_;
        }

        ~FSGenPairFunction() override = default;

        double evaluate(const double distance) override {
            if (distance >= cutoff_) return 0;
            return prefactor_ * pow(1 - distance * cutoff_inverse_, degree_);
        };

        double differentiate(const double distance) override {
            if (distance >= cutoff_) return 0;
            return - prefactor_ * pow(1 - distance * cutoff_inverse_, degree_ - 1) * degree_ * cutoff_inverse_;
        }

    private:
        double cutoff_inverse_;
        double degree_;
    };

    inline std::shared_ptr<FSGenPairFunction> FSGenPairFunction::fromDataNode(const DataNode& params) {
        return std::make_shared<FSGenPairFunction>(
            require(params, "cutoff"),
            require(params, "degree"),
            params.getOrDefault("prefactor", 1.0)
        );
    }

    inline DataNode FSGenPairFunction::serialize() {
        return DataNode{
            {"prefactor", prefactor_},
            {"cutoff", cutoff_},
            {"degree", degree_}
        };
    }
}

#endif //FSGENPAIRFUNCTION_HPP
