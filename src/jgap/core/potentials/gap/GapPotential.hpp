#ifndef JGAP_GAPPOTENTIAL_HPP
#define JGAP_GAPPOTENTIAL_HPP
#include "../Potential.hpp"
#include "component/GapComponent.hpp"

namespace jgap {
    class GapPotential : public Potential {
    public:
        ValuePtr<Potential> optional_external_potential{};
        std::vector<ValuePtr<GapComponent>> components{};

        GapPotential() = default;

        template<typename GapComponentT>
            requires(std::convertible_to<GapComponentT, ValuePtr<GapComponent>>)
        GapPotential(std::initializer_list<GapComponentT> components, ValuePtr<Potential> external = nullptr,
                     const std::vector<Real>& coefficients = {}) :
            optional_external_potential(std::move(external)), components(components.begin(), components.end()) {
            if (!coefficients.empty()) {
                setCoefficients(coefficients);
            }
        }

        template<typename GapComponentT>
        requires (std::convertible_to<GapComponentT, ValuePtr<GapComponent>>)
        void addComponents(std::vector<GapComponentT> new_components) {
            for (auto& comp : new_components) {
                components.push_back(std::move(comp));
            }
        }

        void addComponent(ValuePtr<GapComponent> component);

        void setCoefficients(const std::vector<Real>& new_coefficients);

        AtomicQuantity calculateEnergy(const Atoms& atoms) const override;

        Cutoffs getCutoffs() const override;

        const std::vector<ValuePtr<GapComponent>>& getComponents() const;

        void fillTables(TabulationData& table) const override;

        GapPotential* clone() const override { return new GapPotential(*this); }
    };
}

#endif
