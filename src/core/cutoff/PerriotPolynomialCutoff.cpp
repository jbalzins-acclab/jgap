//
// Created by Jegors Balzins on 9.7.2025.
//
#include "core/cutoff/PerriotPolynomialCutoff.hpp"

namespace jgap {
    PerriotPolynomialCutoff::PerriotPolynomialCutoff(const nlohmann::json &params) {
        _cutoff = params["cutoff"].get<double>();
        if (params.contains("r_min")) {
            _rMin = params["r_min"].get<double>();
        } else {
            _rMin = _cutoff - params["cutoff_transition_width"].get<double>();
        }
        _cutoffWidthInverse = 1.0 / (_cutoff - _rMin);
    }

    nlohmann::json PerriotPolynomialCutoff::serialize() {
        return {
            {"r_min", _rMin},
            {"cutoff", _cutoff}
        };
    }

    PerriotPolynomialCutoff::PerriotPolynomialCutoff(const double rMin, const double cutoff) {
        _rMin = rMin;
        _cutoff = cutoff;
        _cutoffWidthInverse = 1.0 / (_cutoff - _rMin);
    }

    double PerriotPolynomialCutoff::evaluate(double r) {
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
        if (r <= _rMin) return 1.0;
        if (r >= _cutoff) return 0.0;

        const double chi = (r - _rMin) * _cutoffWidthInverse;
        return 1.0 - chi * chi * chi * (6.0 * chi * chi - 15.0 * chi + 10.0);
    }

    double PerriotPolynomialCutoff::differentiate(double r) {
        if (r <= _rMin || r >= _cutoff) return 0.0;

        const double chi = (r - _rMin) * _cutoffWidthInverse;
        const double dchi_dr = _cutoffWidthInverse;

        return dchi_dr * (
                   - 3.0 * chi * chi * (6.0 * chi * chi - 15.0 * chi + 10.0)
                   - chi * chi * chi * (12.0 * chi - 15.0)
               );
    }
}
