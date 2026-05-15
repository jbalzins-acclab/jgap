#ifndef JGAP_REGULARIZATIONRULES_HPP
#define JGAP_REGULARIZATIONRULES_HPP

#include "Regularization.hpp"
#include "core/atomic/neighbours/NeighbourList.hpp"

namespace jgap {
    class RegularizationRules {
    public:
        virtual ~RegularizationRules() = default;
        virtual void fillSigmas(Regularization& sigmas, const Atoms& atoms, const NeighbourList& neighbour_list) = 0;

        /*
        static void fillSigmas(Atoms& atoms, NeighbourList ne,
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
        }*/
    };
}

#endif
