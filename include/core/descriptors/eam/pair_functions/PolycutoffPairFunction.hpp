//
// Created by Jegors Balzins on 19.6.2025.
//

#ifndef POLYCUTOFFPAIRFUNCTION_HPP
#define POLYCUTOFFPAIRFUNCTION_HPP

#include "EamPairFunction.hpp"

namespace jgap {
    class PolycutoffPairFunction : public EamPairFunction {
    public:
        PolycutoffPairFunction(const double cutoff, const double rmin)
            : _cutoff(cutoff), _rmin(rmin) {
            _intervalInverse = 1 / (_cutoff - _rmin);
        }

        ~PolycutoffPairFunction() override = default;

        double evaluate(const double distance) override {
            if (distance >= _cutoff) return 0.0;
            if (distance <= _rmin) return 1.0;

            const double chi = (distance - _rmin) * _intervalInverse;
            return 1.0 - chi * chi * chi * (6 * chi * chi - 15 * chi + 10);
        }

        double differentiate(const double distance) override {
            if (distance >= _cutoff || distance <= _rmin) return 0;

            const double chi = (distance - _rmin) * _intervalInverse;
            const double dchi_dr = _intervalInverse;

            return dchi_dr * chi * chi * ( -30 * chi * chi + 60 * chi - 30);
        }

    private:
        double _cutoff;
        double _rmin;
        double _intervalInverse;
    };
}


#endif //POLYCUTOFFPAIRFUNCTION_HPP
