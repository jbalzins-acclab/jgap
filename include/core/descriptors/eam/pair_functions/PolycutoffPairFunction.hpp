#ifndef POLYCUTOFFPAIRFUNCTION_HPP
#define POLYCUTOFFPAIRFUNCTION_HPP

#include "EamPairFunction.hpp"
#include "io/parse/ParserRegistry.hpp"

namespace jgap {
    class PolycutoffPairFunction : public EamPairFunction, Serializable {
    public:

        SETUP_PARSER_AND_SERIALIZATION(EamPairFunction, PolycutoffPairFunction, polycutoff)

        PolycutoffPairFunction(const double cutoff, const double rMin, const double prefactor = 1.0)
            : EamPairFunction(cutoff, prefactor), _rMin(rMin)
        {
            _intervalInverse = 1.0 / (cutoff_ - _rMin);
        }

        ~PolycutoffPairFunction() override = default;

        double evaluate(const double distance) override {
            if (distance >= cutoff_) return 0.0;
            if (distance <= _rMin) return prefactor_;

            const double chi = (distance - _rMin) * _intervalInverse;
            return prefactor_ * (1.0 - chi * chi * chi * (6 * chi * chi - 15 * chi + 10));
        }

        double differentiate(const double distance) override {
            if (distance >= cutoff_ || distance <= _rMin) return 0;

            const double chi = (distance - _rMin) * _intervalInverse;
            const double dchi_dr = _intervalInverse;

            return prefactor_ * (dchi_dr * chi * chi * ( -30 * chi * chi + 60 * chi - 30));
        }

    private:
        double _rMin;
        double _intervalInverse;
    };

    inline std::shared_ptr<PolycutoffPairFunction> PolycutoffPairFunction::fromDataNode(const DataNode &params) {
        return std::make_shared<PolycutoffPairFunction>(
            require(params, "cutoff"),
            params.getOrDefault(
                "r_min",
            require(params, "cutoff").asDouble() - require(params, "cutoff_transition_width").asDouble()
            ),
            params.getOrDefault("prefactor", 1.0)
        );
    }

    inline DataNode PolycutoffPairFunction::serialize() {
        return DataNode{
            {"prefactor", prefactor_},
            {"cutoff", cutoff_},
            {"r_min", _rMin}
        };
    }
}


#endif //POLYCUTOFFPAIRFUNCTION_HPP
