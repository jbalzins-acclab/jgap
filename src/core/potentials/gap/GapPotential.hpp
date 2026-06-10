#ifndef JGAP_GAPPOTENTIAL_HPP
#define JGAP_GAPPOTENTIAL_HPP
#include "component/GapComponent.hpp"
#include "../Potential.hpp"

namespace jgap {
    class GapPotential : public Potential {
    public:
        std::unique_ptr<Potential> optional_external_potential{};
        std::vector<GapComponent::Ptr> components{};
        //GapPotential() = default; //deserialize

        GapPotential(std::vector<GapComponent::Ptr> components,
                     std::unique_ptr<Potential> external = nullptr,
                     const std::vector<Real> &coefficients = {});

        void setCoefficients(const std::vector<Real>& new_coefficients);

        AtomicQuantity calculateEnergy(const Atoms &atoms) override;

        Cutoffs getCutoffs() override;

        const std::vector<GapComponent::Ptr>& getComponents() const;
    };
}

#endif