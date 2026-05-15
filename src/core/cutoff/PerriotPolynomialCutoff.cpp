//
// Created by Jegors Balzins on 9.7.2025.
//
#include "core/cutoff/PerriotPolynomialCutoff.hpp"

#include "core/Real.hpp"

namespace jgap {

    PerriotPolynomialCutoff::PerriotPolynomialCutoff(const Real r_min, const Real cutoff)
        : r_min(r_min), cutoff(cutoff)
    {
        cutoff_width_inverse = static_cast<Real>(1.0) / (cutoff - r_min);
    }

    /*PerriotPolynomialCutoff::PerriotPolynomialCutoff(const DataNode &params) {
        _cutoff = params["cutoff"].get<Real>();
        if (params.contains("r_min")) {
            _rMin = params["r_min"].get<Real>();
        } else {
            _rMin = _cutoff - params["cutoff_transition_width"].get<Real>();
        }
        _cutoffWidthInverse = 1.0 / (_cutoff - _rMin);
    }

    DataNode PerriotPolynomialCutoff::serialize() {
        return {
            {"r_min", _rMin},
            {"cutoff", _cutoff}
        };
    }*/

    Real PerriotPolynomialCutoff::evaluate(Real r) const {
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
        if (r <= r_min) return 1.0;
        if (r >= cutoff) return 0.0;

        const Real chi = (r - r_min) * cutoff_width_inverse;
        return 1.0 - chi * chi * chi * (6.0 * chi * chi - 15.0 * chi + 10.0);
    }

    Real PerriotPolynomialCutoff::differentiate(Real r) const {
        if (r <= r_min || r >= cutoff) return 0.0;

        const Real chi = (r - r_min) * cutoff_width_inverse;
        const Real dchi_dr = cutoff_width_inverse;

        return dchi_dr * (
                   - 3.0 * chi * chi * (6.0 * chi * chi - 15.0 * chi + 10.0)
                   - chi * chi * chi * (12.0 * chi - 15.0)
               );
    }
}
