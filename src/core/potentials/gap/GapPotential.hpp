#ifndef JGAP_GAPPOTENTIAL_HPP
#define JGAP_GAPPOTENTIAL_HPP
#include "component/GapComponent.hpp"
#include "../Potential.hpp"

namespace jgap {
    class GapPotential : public Potential {
    public:
        ValuePtr<Potential> optional_external_potential{};
        std::vector<ValuePtr<GapComponent>> components{};

        template<typename GapComponentPtr>
        requires (std::convertible_to<GapComponentPtr, ValuePtr<GapComponent>>)
        GapPotential(std::vector<GapComponentPtr> components,
                     ValuePtr<Potential> external = nullptr,
                     const std::vector<Real> &coefficients = {});

        void setCoefficients(const std::vector<Real>& new_coefficients);

        AtomicQuantity calculateEnergy(const Atoms &atoms) const override;

        Cutoffs getCutoffs() const override;

        const std::vector<ValuePtr<GapComponent>>& getComponents() const;

        void tabulate(TabulationData &table) const override;

        std::unique_ptr<Potential> clone() const override {
            return std::make_unique<GapPotential>(*this);
        }
    };

    template<typename GapComponentPtr> requires (std::convertible_to<GapComponentPtr, ValuePtr<GapComponent>>)
    GapPotential::GapPotential(std::vector<GapComponentPtr> components,
                               ValuePtr<Potential> external,
                               const std::vector<Real> &coefficients)
        : optional_external_potential(std::move(external)),
          components(std::move(components)) {

        if (coefficients.empty()) return;
        setCoefficients(coefficients);
    }
}

#endif