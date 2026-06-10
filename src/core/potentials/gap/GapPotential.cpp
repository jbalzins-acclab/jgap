#include "GapPotential.hpp"

namespace jgap {
    GapPotential::GapPotential(std::vector<GapComponent::Ptr> components,
                               std::unique_ptr<Potential> external,
                               const std::vector<Real> &coefficients)
        : optional_external_potential(std::move(external)), components(std::move(components)) {

        if (coefficients.empty()) return;
        setCoefficients(coefficients);
    }

    void GapPotential::setCoefficients(const std::vector<Real> &new_coefficients) {
        if (new_coefficients.empty()) {
            JGAP_LOG_AND_THROW("Cannot set empty coefficients");
        }

        size_t n = std::transform_reduce(
            components.begin(),
            components.end(),
            size_t{0},
            std::plus(),
            [](const auto& ptr) { return ptr->nSparsePoints(); }
        );

        if (n != new_coefficients.size()) {
            JGAP_LOG_AND_THROW("{} sparse points, but {} coefficients", n, new_coefficients.size());
        }

        auto iter = new_coefficients.begin();
        for (auto& component: components) {
            component->setCoefficients(iter);
        }
    }

    AtomicQuantity GapPotential::calculateEnergy(const Atoms &atoms) {
        AtomicQuantity result(atoms.nAtoms());

        const NeighbourList neighbour_list(atoms, getCutoffs().maxOverall());

        if (optional_external_potential != nullptr) {
            result += optional_external_potential->calculateEnergy(atoms);
        }

        for (const auto& component: components) {
            result += component->energy(neighbour_list);
        }

        return result;
    }

    Cutoffs GapPotential::getCutoffs() {
        Cutoffs result{};

        if (optional_external_potential != nullptr) {
            result = optional_external_potential->getCutoffs();
        }

        for (const auto& component: components) {
            result += component->getCutoffs();
        }

        return result;
    }

    const std::vector<GapComponent::Ptr>& GapPotential::getComponents() const {
        return components;
    }

    /*
    DataNode GapPotential::serializeWithoutType() {
        DataNode result{
            {"coefficients", coefficients}
        };

        if (optional_external_potential != nullptr) {
            result["external_potential"] = optional_external_potential->serialize();
        }

        result["components"] = DataNode::array();
        for (size_t i = 0; i < components.size(); i++) {
            result["components"][i] = components[i]->serialize();
        }

        return result;
    }

    void GapPotential::deserializeNoTypeCheck(const DataNode &serialized) {
        auto& unparsed_coeff = REQUIRE(serialized, "coefficients").asArray();
        coefficients.resize(unparsed_coeff.size());

        std::ranges::transform(unparsed_coeff, coefficients.begin(),
                               [](const DataNode& v) -> Real {return v.asDouble();});

        if (serialized.contains("external_potential")) {
            if (optional_external_potential == nullptr) {
                JGAP_LOG_AND_THROW("Found external potential in serialized data, "
                                   "but optional_external_potential was not specified");
            }
            optional_external_potential->deserialize(serialized["external_potential"]);
        } else {
            if (optional_external_potential != nullptr) {
                JGAP_LOG_AND_THROW("Found no external potential in serialized data, "
                                   "but optional_external_potential of type={} was specified",
                                   optional_external_potential->getTypeId());
            }
        }
    }
    */
}
