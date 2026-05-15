#ifndef JGAP_GAPPOTENTIAL_HPP
#define JGAP_GAPPOTENTIAL_HPP
#include "GapComponent.hpp"
#include "../Potential.hpp"

namespace jgap {
    class GapPotential : public Potential {
    public:
        std::shared_ptr<Potential> optional_external_potential{};
        std::vector<IGapComponent::Ptr> components{};
        //GapPotential() = default; //deserialize

        GapPotential(const std::vector<IGapComponent::Ptr> &components,
                     const std::shared_ptr<Potential> &external = nullptr,
                     const std::vector<Real> &coefficients = {});

        void setCoefficients(const std::vector<Real>& new_coefficients);

        AtomicQuantity energy(const Atoms &atoms) override;

        Cutoffs getCutoffs() override;

        std::string getTypeId() {
            return "GapPotential";
        }

        /*
    protected:
        DataNode serializeWithoutType() override;
        void deserializeNoTypeCheck(const DataNode &serialized) override;
*/
    private:
        std::vector<Real> coefficients{};
    };
}

#endif