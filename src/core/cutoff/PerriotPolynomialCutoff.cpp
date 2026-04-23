//
// Created by Jegors Balzins on 9.7.2025.
//
#include "core/cutoff/PerriotPolynomialCutoff.hpp"

namespace jgap {

    PerriotPolynomialCutoff::PerriotPolynomialCutoff(const double rMin, const double cutoff) {
        r_min_ = rMin;
        cutoff_ = cutoff;
        cutoff_width_inverse_ = 1.0 / (cutoff_ - r_min_);
    }

    /*PerriotPolynomialCutoff::PerriotPolynomialCutoff(const DataNode &params) {
        _cutoff = params["cutoff"].get<double>();
        if (params.contains("r_min")) {
            _rMin = params["r_min"].get<double>();
        } else {
            _rMin = _cutoff - params["cutoff_transition_width"].get<double>();
        }
        _cutoffWidthInverse = 1.0 / (_cutoff - _rMin);
    }

    DataNode PerriotPolynomialCutoff::serialize() {
        return {
            {"r_min", _rMin},
            {"cutoff", _cutoff}
        };
    }*/

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
        if (r <= r_min_) return 1.0;
        if (r >= cutoff_) return 0.0;

        const double chi = (r - r_min_) * cutoff_width_inverse_;
        return 1.0 - chi * chi * chi * (6.0 * chi * chi - 15.0 * chi + 10.0);
    }

    double PerriotPolynomialCutoff::differentiate(double r) {
        if (r <= r_min_ || r >= cutoff_) return 0.0;

        const double chi = (r - r_min_) * cutoff_width_inverse_;
        const double dchi_dr = cutoff_width_inverse_;

        return dchi_dr * (
                   - 3.0 * chi * chi * (6.0 * chi * chi - 15.0 * chi + 10.0)
                   - chi * chi * chi * (12.0 * chi - 15.0)
               );
    }
}
