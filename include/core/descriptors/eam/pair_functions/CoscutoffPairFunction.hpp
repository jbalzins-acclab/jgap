#ifndef COSCUTOFFPAIRFUNCTION_HPP
#define COSCUTOFFPAIRFUNCTION_HPP

#include <math.h>

#include "EamPairFunction.hpp"
#include "io/parse/ParserRegistry.hpp"

namespace jgap {
    class CoscutoffPairFunction : public EamPairFunction {
    public:
        CoscutoffPairFunction(const nlohmann::json& params) {
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

        CoscutoffPairFunction(const double cutoff, const double rmin)
            : _cutoff(cutoff), _rmin(rmin) {
            _intervalInverse = 1.0 / (_cutoff - _rmin);
        }

        ~CoscutoffPairFunction() override = default;

        double evaluate(const double distance) override {
            if (distance >= _cutoff) return 0.0;
            if (distance <= _rmin) return 1.0;

            const double chi = (distance - _rmin) * _intervalInverse;
            return _prefactor * 0.5 * (1 + cos(M_PI * chi));
        };

        double differentiate(const double distance) override {
            if (distance >= _cutoff || distance <= _rmin) return 0;

            const double chi = (distance - _rmin) * _intervalInverse;
            const double dchi_dr = _intervalInverse;

            return -_prefactor * dchi_dr * 0.5 * M_PI * sin(M_PI * chi);
        }

        string getType() override { return "coscutoff"; }

    private:
        double _cutoff;
        double _rmin;
        double _intervalInverse;
    };

    REGISTER_PARSER("coscutoff", EamPairFunction, CoscutoffPairFunction)
}

#endif //COSCUTOFFPAIRFUNCTION_HPP
