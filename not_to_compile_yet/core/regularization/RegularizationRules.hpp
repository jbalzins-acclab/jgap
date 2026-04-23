#ifndef JGAP_REGULARIZATIONRULES_HPP
#define JGAP_REGULARIZATIONRULES_HPP

#include "data/atomic/AtomicStructure.hpp"

namespace jgap {
    class RegularizationRules {
    public:
        struct Conventions {
            static constexpr std::string ENERGY_PER_ATOM = "energy_per_atom";
            static constexpr std::string FORCE_COMPONENT = "force_component";
            static constexpr std::string FORCE_MAGNITUDE = "force_magnitude";
            static constexpr std::string VIRIAL_PER_ATOM = "virial_per_atom";
            static constexpr std::string VIRIAL_ISO_PER_ATOM = "virial_iso_per_atom";
            static constexpr std::string VIRIAL_ANISO_PER_ATOM = "virial_aniso_per_atom";
        };

        virtual ~RegularizationRules() = default;
        virtual void fillSigmas(std::vector<AtomicStructure>& structures) = 0;

    protected:
        static void fillSigmas(AtomicStructure& structure, double energyPerAtom,
                                const Vector3 &force, double virialIso, double virialAniso) {

            if (!structure.energySigmaInverse.has_value()) {
                structure.energySigmaInverse = 1.0 / (energyPerAtom * pow(structure.size(), 0.5));
            }

            if (!structure.forceSigmasInverse.has_value()) {
                structure.forceSigmasInverse = std::vector(
                    structure.size(),  Vector3{1.0 / force.x, 1.0 / force.y, 1.0 / force.z}
                    );
            }

            const double dV_iso = 1.0 / (virialIso * pow(structure.size(), 0.5));
            const double dV_aniso = 1.0 / (virialIso * pow(structure.size(), 0.5));
            if (!structure.virialSigmasInverse.has_value()) {
                structure.virialSigmasInverse = {
                    Vector3{dV_iso, dV_aniso, dV_aniso},
                    Vector3{dV_aniso, dV_iso, dV_aniso},
                    Vector3{dV_aniso, dV_aniso, dV_iso}
                };
            }
        }
    };
}

#endif
