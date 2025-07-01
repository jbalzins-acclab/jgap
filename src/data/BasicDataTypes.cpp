#include "data/BasicDataTypes.hpp"

namespace jgap {
    void AtomicStructure::set(const PotentialPrediction &prediction) {
        adjust(prediction, false, true);
    }

    void AtomicStructure::adjust(const PotentialPrediction &prediction,
                                 const bool subtract,
                                 const bool setEmpty) {

        const double sign = subtract ? -1.0 : 1.0;

        if (prediction.energy.has_value() && (this->energy.has_value() || setEmpty)) {
            energy = energy.value_or(0) + sign * prediction.energy.value();
        }

        if (prediction.forces.has_value() && (setEmpty || atoms[0].force.has_value())) {
            if (this->atoms.size() != prediction.forces.value().size()) {
                const string errMsg = format(
                    "Found force {} predictions for a {} atom system",
                    prediction.forces.value().size(),
                    this->atoms.size()
                );
                if (Logger::logger != nullptr) {
                    Logger::logger->error(errMsg);
                }
                throw runtime_error(errMsg);
            }

            for (size_t i = 0; i < this->atoms.size(); i++) {
                atoms[i].force = atoms[i].force.value_or(Vector3{0,0,0})
                                 + prediction.forces.value()[i] * sign;
            }
        }

        if (prediction.virials.has_value() && (this->virials.has_value() || setEmpty)) {
            if (!this->virials.has_value()) {
                this->virials = array{array{0.0,0.0,0.0}, array{0.0,0.0,0.0}, array{0.0,0.0,0.0}};
            }
            for (size_t i = 0; i < 3; i++) {
                for (size_t j = 0; j < 3; j++) {
                    (*this->virials)[i][j] += sign * prediction.virials.value()[i][j];
                }
            }
        }
    }

    AtomicStructure AtomicStructure::repeat(size_t a, size_t b, size_t c) {
        auto cpy = AtomicStructure(*this);

        cpy.latticeVectors[0] = cpy.latticeVectors[0] * a;
        cpy.latticeVectors[1] = cpy.latticeVectors[1] * b;
        cpy.latticeVectors[2] = cpy.latticeVectors[2] * c;

        cpy.atoms = {};
        for (size_t i = 0; i < a; i++) {
            for (size_t j = 0; j < b; j++) {
                for (size_t k = 0; k < c; k++) {
                    for (auto& atom: this->atoms) {
                        auto newAtom = AtomData(atom);
                        newAtom.position = newAtom.position
                            + this->latticeVectors[0] * i
                            + this->latticeVectors[1] * j
                            + this->latticeVectors[2] * k;
                        cpy.atoms.push_back(newAtom);
                    }
                }
            }
        }

        return cpy;
    }
}
