//
// Created by Jegors Balzins on 19.6.2025.
//

#ifndef POLYCUTOFFPAIRFUNCTION_HPP
#define POLYCUTOFFPAIRFUNCTION_HPP

#include "EamPairFunction.hpp"
#include "io/parse/ParserRegistry.hpp"

namespace jgap {
    class PolycutoffPairFunction : public EamPairFunction {
    public:
        explicit PolycutoffPairFunction(nlohmann::json params) {
            _cutoff = params["cutoff"];
            if (params.contains("r_min")) {
                _rmin = params["r_min"];
            } else if (params.contains("cutoff_transition_width")) {
                _rmin = _cutoff - params["cutoff_transition_width"].get<double>();
            } else {
                _rmin = 0;
            }
            _prefactor = params.value("prefactor", 1.0);
            _intervalInverse = 1.0 / (_cutoff - _rmin);
        }
        nlohmann::json serialize() override {
            return {
                {"prefactor", _prefactor},
                {"cutoff", _cutoff},
                {"r_min", _rmin}
            };
        }

        PolycutoffPairFunction(const double cutoff, const double rmin)
            : _cutoff(cutoff), _rmin(rmin) {
            _intervalInverse = 1.0 / (_cutoff - _rmin);
        }

        ~PolycutoffPairFunction() override = default;

        double evaluate(const double distance) override {
            if (distance >= _cutoff) return 0.0;
            if (distance <= _rmin) return 1.0;

            const double chi = (distance - _rmin) * _intervalInverse;
            return _prefactor * (1.0 - chi * chi * chi * (6 * chi * chi - 15 * chi + 10));
        }

        double differentiate(const double distance) override {
            if (distance >= _cutoff || distance <= _rmin) return 0;

            const double chi = (distance - _rmin) * _intervalInverse;
            const double dchi_dr = _intervalInverse;

            return _prefactor * (dchi_dr * chi * chi * ( -30 * chi * chi + 60 * chi - 30));
        }

        string getType() override { return "polycutoff"; }

    private:
        double _cutoff;
        double _rmin;
        double _intervalInverse;
    };

    REGISTER_PARSER("polycutoff", EamPairFunction, PolycutoffPairFunction)
}


#endif //POLYCUTOFFPAIRFUNCTION_HPP
