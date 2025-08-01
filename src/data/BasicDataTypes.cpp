#include "data/BasicDataTypes.hpp"

#include <format>

namespace jgap {
    void AtomicStructure::setEnergyData(const PotentialPrediction &prediction) {
        energy = {};
        forces = {};
        virials = {};
        adjust(prediction, false, true);
    }

    void AtomicStructure::adjust(const PotentialPrediction &prediction,
                                 const bool subtract,
                                 const bool setEmpty) {

        const double sign = subtract ? -1.0 : 1.0;

        if (prediction.energy.has_value() && (this->energy.has_value() || setEmpty)) {
            energy = energy.value_or(0) + sign * prediction.energy.value();
        }

        if (prediction.forces.has_value() && (setEmpty || forces.has_value())) {
            if (size() != prediction.forces.value().size()) {
                const string errMsg = format(
                    "Found force {} predictions for a {} atom system",
                    prediction.forces.value().size(),
                    size()
                );
                if (CurrentLogger::get() != nullptr) {
                    CurrentLogger::get()->error(errMsg);
                }
                throw runtime_error(errMsg);
            }

            if (!forces.has_value()) {
                forces = vector(size(), Vector3{0.0, 0.0, 0.0});
            }
            for (size_t i = 0; i < size(); i++) {
                (*forces)[i] += prediction.forces.value()[i] * sign;
            }
        }

        if (prediction.virials.has_value() && (this->virials.has_value() || setEmpty)) {
            if (!this->virials.has_value()) {
                this->virials = array{Vector3{0.0,0.0,0.0}, Vector3{0.0,0.0,0.0}, Vector3{0.0,0.0,0.0}};
            }
            for (size_t i = 0; i < 3; i++) {
                (*this->virials)[i] = (*this->virials)[i] + prediction.virials.value()[i] * sign;
            }
        }
    }

    AtomicStructure AtomicStructure::repeat(size_t a, size_t b, size_t c) {
        auto cpy = AtomicStructure{};
        cpy.configType = this->configType;

        cpy.energy = this->energy.transform([&](const double val) -> double {
            return val * static_cast<double>(a) * static_cast<double>(b) * static_cast<double>(c);
        });

        cpy.lattice[0] = this->lattice[0] * static_cast<double>(a);
        cpy.lattice[1] = this->lattice[1] * static_cast<double>(b);
        cpy.lattice[2] = this->lattice[2] * static_cast<double>(c);

        cpy.species = {};
        cpy.positions = {};
        cpy.forces = {};
        cpy.forceSigmasInverse = {};

        if (forces.has_value()) {
            cpy.forces = vector<Vector3>{};
        }
        for (size_t i = 0; i < a; i++) {
            for (size_t j = 0; j < b; j++) {
                for (size_t k = 0; k < c; k++) {
                    for (const auto& atom: *this) {
                        cpy.species.push_back(atom.species());
                        cpy.positions.push_back(
                            atom.position()
                                + this->lattice[0] * static_cast<double>(i)
                                + this->lattice[1] * static_cast<double>(j)
                                + this->lattice[2] * static_cast<double>(k)
                            );
                        if (forces.has_value()) {
                            cpy.forces->push_back(atom.force());
                        }
                    }
                }
            }
        }

        return cpy;
    }

    double AtomicStructure::volume() const {
        return abs(lattice[0].cross(lattice[1]).dot(lattice[2]));
    }
}
