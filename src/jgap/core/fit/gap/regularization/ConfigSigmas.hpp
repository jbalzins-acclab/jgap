#ifndef JGAP_CONFIGSIGMAS_HPP
#define JGAP_CONFIGSIGMAS_HPP

#include "jgap/core/Vector3.hpp"
#include "jgap/core/atomic/energy/Virials.hpp"

namespace jgap {
    struct ConfigSigmas {
        Real energy;
        Vector3 force;
        Virials virials;

        ConfigSigmas(Real e) :
            energy(e),
            force{e * 50.0_r, e * 50.0_r, e * 50.0_r},
            virials{e * 100.0_r, e * 100.0_r, e * 100.0_r, e * 100.0_r, e * 100.0_r, e * 100.0_r} {}

        ConfigSigmas(Real e, Real f, Real v) : energy(e), force{f, f, f}, virials{v, v, v, v, v, v} {}

        ConfigSigmas(Real e, Real f, Real v_iso, Real v_aniso) :
            energy(e), force{f, f, f}, virials{v_iso, v_aniso, v_aniso, v_iso, v_aniso, v_iso} {}

        ConfigSigmas(Real e, Vector3 f, Virials v) : energy(e), force(f), virials(v) {}
    };
}

#endif
