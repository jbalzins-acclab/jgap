//
// Created by Jegors Balzins on 19.6.2025.
//

#ifndef COSCUTOFFPAIRFUNCTION_HPP
#define COSCUTOFFPAIRFUNCTION_HPP

#include <math.h>

#include "EamPairFunction.hpp"

namespace jgap {
    class CoscutoffPairFunction : public EamPairFunction{
    public:
        CoscutoffPairFunction(const double cutoff, const double rmin)
            : _cutoff(cutoff), _rmin(rmin) {
            _intervalInverse = 1 / (_cutoff - _rmin);
        }

        ~CoscutoffPairFunction() override = default;

        double evaluate(const double distance) override {
            if (distance >= _cutoff) return 0.0;
            if (distance <= _rmin) return 1.0;

            const double chi = (distance - _rmin) * _intervalInverse;
            return 0.5 * (1 + cos(M_PI * chi));
        };

        double differentiate(const double distance) override {
            if (distance >= _cutoff || distance <= _rmin) return 0;

            const double chi = (distance - _rmin) * _intervalInverse;
            const double dchi_dr = _intervalInverse;

            return - dchi_dr * 0.5 * M_PI * sin(M_PI * chi);
        }
    private:
        double _cutoff;
        double _rmin;
        double _intervalInverse;
    };
}

#endif //COSCUTOFFPAIRFUNCTION_HPP
