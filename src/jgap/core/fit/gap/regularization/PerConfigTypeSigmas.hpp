#ifndef JGAP_PERCONFIGTYPESIGMAS_HPP
#define JGAP_PERCONFIGTYPESIGMAS_HPP

#include "jgap/core/Vector3.hpp"
#include "jgap/core/atomic/energy/Virials.hpp"

namespace jgap {
    struct PerConfigTypeSigmas {
        double energy;
        Vector3 force;
        Virials virials;

        PerConfigTypeSigmas(double e) :
            energy(e),
            force{e * 50.0, e * 50.0, e * 50.0},
            virials{e * 100.0, e * 100.0, e * 100.0, e * 100.0, e * 100.0, e * 100.0} {}

        PerConfigTypeSigmas(double e, double f, double v) : energy(e), force{f, f, f}, virials{v, v, v, v, v, v} {}

        PerConfigTypeSigmas(double e, double f, double v_iso, double v_aniso) :
            energy(e), force{f, f, f}, virials{v_iso, v_aniso, v_aniso, v_iso, v_aniso, v_iso} {}

        PerConfigTypeSigmas(double e, Vector3 f, Virials v) : energy(e), force(f), virials(v) {}
    };
}

#endif
