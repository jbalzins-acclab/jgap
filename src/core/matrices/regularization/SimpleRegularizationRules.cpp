#include "core/matrices/regularization/SimpleRegularizationRules.hpp"

namespace jgap {

    std::shared_ptr<SimpleRegularizationRules> SimpleRegularizationRules::fromDataNode(const DataNode &params) {
        double e = require(params, Conventions::ENERGY_PER_ATOM);
        double f = params.getOrDefault(Conventions::FORCE_COMPONENT, e * 50.0);

        double vIso, vAniso;
        if (params.contains(Conventions::VIRIAL_ISO_PER_ATOM)
            || params.contains(Conventions::VIRIAL_ANISO_PER_ATOM)) {

            vIso = require(params, Conventions::VIRIAL_ISO_PER_ATOM);
            vAniso = require(params, Conventions::VIRIAL_ANISO_PER_ATOM);

        } else if (params.contains(Conventions::VIRIAL_PER_ATOM)) {
            vIso = vAniso = require(params, Conventions::VIRIAL_PER_ATOM);
        } else {
            vIso = vAniso = e * 100;
        }

        double liquidMultiplier = require(params, "liquid_multiplier");
        double srMultiplier = require(params, "short_range_multiplier");

        return std::make_shared<SimpleRegularizationRules>(e, f, vIso, vAniso, liquidMultiplier, srMultiplier);
    }

    SimpleRegularizationRules::SimpleRegularizationRules(double energyPerAtom, double forceComponent,
                                                         double virialIsoPerAtom, double virialAnisoPerAtom,
                                                         double liquidMultiplier, double shortRangeMultiplier)
        : _defaultEnergyPerAtom(energyPerAtom), _defaultForceComponent(forceComponent),
          _defaultVirialIsoPerAtom(virialIsoPerAtom), _defaultVirialAnisoPerAtom(virialAnisoPerAtom),
          _liquidMultiplier(liquidMultiplier), _shortRangeMultiplier(shortRangeMultiplier) {
    }

    void SimpleRegularizationRules::fillSigmas(std::vector<AtomicStructure> &structures) {
        for (auto& structure : structures) {
            fillSigmas(structure);
        }
    }

    void SimpleRegularizationRules::fillSigmas(AtomicStructure &structure) {
        double multiplier = 1.0;
        std::string ct = "default";
        if (structure.properties.contains("config_type")) {
            ct = structure.properties["config_type"];
        }

        if (ct == "isolated_atom") {
            multiplier = 0.001;
        }

        if (ct.contains("liquid") || ct.contains("melt")) {
            multiplier = _liquidMultiplier;
        }

        if (ct.contains("short") || ct.contains("traj") || ct.contains("low_volume")
            || ct.contains("dimer") || ct.contains("trimer")) {
            multiplier = _shortRangeMultiplier;
        }

        RegularizationRules::fillSigmas(
            structure, multiplier * _defaultEnergyPerAtom,
            Vector3{_defaultForceComponent, _defaultForceComponent, _defaultForceComponent} * multiplier,
            multiplier * _defaultVirialIsoPerAtom, multiplier * _defaultVirialAnisoPerAtom
            );
    }
}
