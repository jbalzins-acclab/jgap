#ifndef PERRIOTPOLYNOMIALCUTOFF_HPP
#define PERRIOTPOLYNOMIALCUTOFF_HPP

#include "CutoffFunction.hpp"

namespace jgap {
    class PerriotPolynomialCutoff : public CutoffFunction {
    public:
        PerriotPolynomialCutoff(const double cutoffLower, const double cutoffUpper) {
            _cutoffLower = cutoffLower;
            _cutoffUpper = cutoffUpper;
            _cutoffWidthInverse = 1.0 / (_cutoffUpper - _cutoffLower);
        }

        double evaluate(double r) override {
            /*
            # Perriot polynomial cutoff
            if r < rmin:
                y = 1.0
            elif r > rmax:
                y = 0.0
            else:
                chi = (r - rmin) / (rmax - rmin)
                y = 1.0 - chi**3 * (6.0*chi**2 - 15*chi + 10.0)
                */
            if (r <= _cutoffLower) return 1.0;
            if (r >= _cutoffUpper) return 0.0;

            const double chi = (r - _cutoffLower) * _cutoffWidthInverse;
            return 1.0 - chi * chi * chi * (6 * chi * chi - 15 * chi + 10);
        }

        double differentiate(double r) override {
            if (r <= _cutoffLower || r >= _cutoffUpper) return 0.0;

            const double chi = (r - _cutoffLower) * _cutoffWidthInverse;
            const double dchi_dr = _cutoffWidthInverse;

            return dchi_dr * (
                    - 3 * chi * chi * (6 * chi * chi - 15 * chi + 10)
                    - chi * chi * chi * (12 * chi - 15)
                    );
        }

    private:
        double _cutoffUpper;
        double _cutoffLower;
        double _cutoffWidthInverse;
    };
}

#endif //PERRIOTPOLYNOMIALCUTOFF_HPP
