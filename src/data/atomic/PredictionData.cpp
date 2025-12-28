#include "../../../include/data/atomic/PredictionData.hpp"

namespace jgap {
    Predictions Predictions::operator+(const Predictions &other) const {

        std::optional<double> _energy;
        if (energy.has_value() || other.energy.has_value()) {
            _energy = energy.value_or(0.0) + other.energy.value_or(0.0);
        }

        std::optional<std::vector<Vector3>> _forces;
        if (forces.has_value() && other.forces.has_value()) {
            _forces = std::vector<Vector3>(forces.value().size());
            for (size_t i = 0; i < forces.value().size(); i++) {
                _forces->at(i) = forces.value()[i] + other.forces.value()[i];
            }
        } else if (forces.has_value()) {
            _forces = forces.value();
        } else if (other.forces.has_value()) {
            _forces = other.forces.value();
        }

        std::optional<std::array<Vector3, 3>> _virials;
        if (virials.has_value() && other.virials.has_value()) {
            _virials = virials;
            (*_virials)[0] += other.virials.value()[0];
            (*_virials)[1] += other.virials.value()[1];
            (*_virials)[2] += other.virials.value()[2];
        } else if (virials.has_value()) {
            _virials = virials.value();
        } else if (other.virials.has_value()) {
            _virials = other.virials.value();
        }

        return Predictions{
            _energy, _forces, _virials
        };
    }
}