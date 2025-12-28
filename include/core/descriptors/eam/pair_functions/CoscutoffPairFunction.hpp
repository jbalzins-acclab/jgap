#ifndef COSCUTOFFPAIRFUNCTION_HPP
#define COSCUTOFFPAIRFUNCTION_HPP

#include <cmath>

#include "EamPairFunction.hpp"
#include "io/parse/ParserRegistry.hpp"

namespace jgap {
    class CoscutoffPairFunction : public EamPairFunction, Serializable {
    public:
        SETUP_PARSER_AND_SERIALIZATION(EamPairFunction, CoscutoffPairFunction, coscutoff)

        ~CoscutoffPairFunction() override = default;
        CoscutoffPairFunction(const double cutoff, const double rMin, const double prefactor = 1.0) : r_min_(rMin) {
            cutoff_ = cutoff;
            interval_inverse_ = 1.0 / (cutoff_ - r_min_);
            prefactor_ = prefactor;
        }

        double evaluate(const double distance) override {
            if (distance >= cutoff_) return 0.0;
            if (distance <= r_min_) return prefactor_;

            const double chi = (distance - r_min_) * interval_inverse_;
            return prefactor_ * 0.5 * (1 + cos(M_PI * chi));
        };

        double differentiate(const double distance) override {
            if (distance >= cutoff_ || distance <= r_min_) return 0;

            const double chi = (distance - r_min_) * interval_inverse_;
            const double dchi_dr = interval_inverse_;

            return -prefactor_ * dchi_dr * 0.5 * M_PI * sin(M_PI * chi);
        }

    private:
        double r_min_;
        double interval_inverse_;
    };

    inline DataNode CoscutoffPairFunction::serialize() {
        return {
            { "prefactor", prefactor_},
            { "cutoff_", cutoff_},
            { "r_min", r_min_},
        };
    }

    inline std::shared_ptr<CoscutoffPairFunction> CoscutoffPairFunction::fromDataNode(const DataNode &params) {

        double r_min = 0.0;
        if (params.type == DataNode::Type::OBJECT) {
            if (params.contains("r_min")) {
                r_min = require(params, "r_min").asDouble();
            } else if (params.contains("cutoff_transition_width")) {
                r_min = require(params, "cutoff").asDouble() - require(params, "cutoff_transition_width").asDouble();
            } else {
                r_min = 0.0;
            }
        }

        return std::make_shared<CoscutoffPairFunction>(
            require(params, "cutoff"),
            r_min,
            params.getOrDefault("prefactor", 1.0)
        );
    }
}

#endif //COSCUTOFFPAIRFUNCTION_HPP
