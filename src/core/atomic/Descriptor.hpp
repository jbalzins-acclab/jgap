#ifndef JGAP_DESCRIPTOR_HPP
#define JGAP_DESCRIPTOR_HPP

#include "geometry/Vector3.hpp"
#include "core/atomic/energy/Virials.hpp"

#include <array>
#include <vector>

#include "core/CalculationType.hpp"

namespace jgap {

    struct GradientData {
        constexpr static size_t UNKNOWN_DEPENDENCIES = 0;

        size_t wrt_atom_index{};
        Vector3 value{};
    };

    // Primary template is not defined, forcing specialization.
    template<
        size_t Dim,
        size_t DependencyAtoms = GradientData::UNKNOWN_DEPENDENCIES,
        CalculationType Type = CalculationType::ValueOnly
    >
    struct Descriptor;

    // Specializations for ValueOnly
    template<size_t Dim>
    struct Descriptor<Dim> {
        std::array<Real, Dim> value;

        Descriptor& operator+=(const Descriptor& other) {
            for (size_t i = 0; i < Dim; i++) {
                value[i] += other.value[i];
            }
            return *this;
        }

        Descriptor operator+(const Descriptor& other) const {
            Descriptor result = *this;
            result += other;
            return result;
        }
    };

    template<size_t Dim, size_t DependencyAtoms>
    struct Descriptor<Dim, DependencyAtoms, CalculationType::ValueOnly> : Descriptor<Dim> {
        // Inherits `value` and operators from Descriptor<Dim>
    };

    // Specialization for a fixed number of gradients
    template<size_t Dim, size_t DependencyAtoms>
    struct Descriptor<Dim, DependencyAtoms, CalculationType::WithGradients> : Descriptor<Dim, DependencyAtoms> {

        using Descriptor<Dim>::value;

        std::array<Virials, Dim> virials{};
        std::array<std::array<GradientData, DependencyAtoms>, Dim> gradients{};

        Descriptor& operator+=(const Descriptor& other) {
            Descriptor<Dim>::operator+=(other); // Call base operator for `value`
            for (size_t i = 0; i < Dim; i++) {
                virials[i] += other.virials[i];
                for (size_t j = 0; j < DependencyAtoms; j++) {
                    gradients[i][j].value += other.gradients[i][j].value;
                    // Assuming wrt_atom_index are matching
                }
            }
            return *this;
        }

        Descriptor operator+(const Descriptor& other) const {
            Descriptor result = *this;
            result += other;
            return result;
        }
    };

    // Specialization for unknown number of gradients
    template<size_t Dim>
    struct Descriptor<Dim, GradientData::UNKNOWN_DEPENDENCIES, CalculationType::WithGradients> : Descriptor<Dim> {

        using Descriptor<Dim>::value;

        std::array<Virials, Dim> virials;
        std::array<std::vector<Vector3>, Dim> gradients;

        explicit Descriptor(size_t box_n_atoms) {
            for (size_t dim = 0; dim < Dim; dim++) {
                gradients[dim].resize(box_n_atoms);
            }
        }
        Descriptor() = default;

        Descriptor& operator+=(const Descriptor& other) {
            Descriptor<Dim>::operator+=(other); // Call base operator for `value`
            for (size_t i = 0; i < Dim; i++) {
                virials[i] += other.virials[i];
                for (size_t j = 0; j < gradients[i].size(); j++) {
                    gradients[i][j] += other.gradients[i][j];
                }
            }
            return *this;
        }

        Descriptor operator+(const Descriptor& other) const {
            Descriptor result = *this;
            result += other;
            return result;
        }
    };
}

#endif
