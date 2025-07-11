//
// Created by Jegors Balzins on 19.6.2025.
//

#ifndef FSGENPAIRFUNCTION_HPP
#define FSGENPAIRFUNCTION_HPP

#include "EamPairFunction.hpp"
#include "io/parse/ParserRegistry.hpp"

namespace jgap {
    class FSGenPairFunction : public EamPairFunction {
    public:
        FSGenPairFunction(const nlohmann::json& params) {
            _cutoff = params["cutoff"];
            _degree = params["degree"];
            _prefactor = params.value("prefactor", 1.0);
            _cutoffInverse = 1.0 / _cutoff;
        }

        nlohmann::json serialize() override {
            return {
                    {"prefactor", _prefactor},
                    {"cutoff", _cutoff},
                    {"degree", _degree}
            };
        }

        FSGenPairFunction(const double cutoff, const double degree) : _cutoff(cutoff), _degree(degree) {
            _cutoffInverse = 1.0 / _cutoff;
        }

        ~FSGenPairFunction() override = default;

        double evaluate(const double distance) override {
            if (distance >= _cutoff) return 0;
            return _prefactor * pow(1 - distance * _cutoffInverse, _degree);
        };

        double differentiate(const double distance) override {
            if (distance >= _cutoff) return 0;
            return - _prefactor * pow(1 - distance * _cutoffInverse, _degree - 1) * _degree * _cutoffInverse;
        }

        string getType() override { return "fsgen"; }

    private:
        double _cutoff;
        double _cutoffInverse;
        double _degree;
    };

    REGISTER_PARSER("fsgen", EamPairFunction, FSGenPairFunction)
}

#endif //FSGENPAIRFUNCTION_HPP
