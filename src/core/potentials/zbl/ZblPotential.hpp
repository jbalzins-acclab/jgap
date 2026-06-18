#ifndef ZBLPOTENTIAL_HPP
#define ZBLPOTENTIAL_HPP

#include <set>
#include <array>
#include <map>
#include <istream>

#include "core/cutoff/CutoffFunction.hpp"
#include "core/cutoff/PerriotPolynomialCutoff.hpp"
#include "core/potentials/Potential.hpp"
#include "../../../serialization/SerializationRegistry.hpp"
#include "core/atomic/Atoms.hpp"
#include "core/atomic/species/SpeciesSet.hpp"
#include "io/log/CurrentLogger.hpp"

namespace jgap {

#if defined(__has_embed)
    enum class EmbeddedZBLCoeffDataset {
        DMOL,
        NLH,
        MP2
    };
#endif

    class ZblPotential : public Potential {
    public:
        static constexpr Real DefaultZblCutoff = 2.2;
        static constexpr Real DefaultZblCutoffTransitionWidth = 1.2;
        static constexpr Real Epsilon0_F_per_m = 8.854187817e-12;
        static constexpr Real ElectronCharge_C = 1.60217657e-19;
        static constexpr Real CoulombConstant_eV_Ang = ElectronCharge_C / (4.0 * M_PI * Epsilon0_F_per_m * 1e-10);

#if defined(__has_embed)
        ZblPotential(const std::set<SpeciesSet<2, Symmetric>>& species,
                     EmbeddedZBLCoeffDataset embedded_dataset = EmbeddedZBLCoeffDataset::DMOL,
                     Real cutoff = DefaultZblCutoff,
                     Real cutoff_transition_width = DefaultZblCutoffTransitionWidth);

        ZblPotential(const std::vector<Atoms>& training_data,
                     EmbeddedZBLCoeffDataset embedded_dataset = EmbeddedZBLCoeffDataset::DMOL,
                     Real cutoff = DefaultZblCutoff,
                     Real cutoff_transition_width = DefaultZblCutoffTransitionWidth);
#endif

        ZblPotential(std::istream& custom_dataset,
                     const std::set<SpeciesSet<2, Symmetric>>& species,
                     Real cutoff = DefaultZblCutoff,
                     Real cutoff_transition_width = DefaultZblCutoffTransitionWidth);

        ZblPotential(std::istream& custom_dataset,
                     const std::vector<Atoms>& training_data,
                     Real cutoff = DefaultZblCutoff,
                     Real cutoff_transition_width = DefaultZblCutoffTransitionWidth);

        // z1_z2 and a_inverse are deduced from the species in each pair.
        ZblPotential(const std::map<SpeciesSet<2, Symmetric>, std::array<Real, 6>>& coefficients,
                     Real cutoff = DefaultZblCutoff,
                     Real cutoff_transition_width = DefaultZblCutoffTransitionWidth);

        std::map<SpeciesSet<2, Symmetric>, std::array<Real, 6>> getCoefficients() const;

        Real getCutoff() const { return cutoff; }
        Real getCutoffTransitionWidth() const { return cutoff_transition_width; }

        AtomicQuantity calculateEnergy(const Atoms &atoms) const override;

        std::array<Real, 2> energyAndDerivative(const SpeciesSet<2, Symmetric>& species_pair, Real r) const;

        Cutoffs getCutoffs() const override { return {{2u, cutoff}}; }

        void fillTables(TabulationData &tables) const override;

        std::unique_ptr<Potential> clone() const override {
            return std::make_unique<ZblPotential>(*this);
        }

    private:
        static constexpr Real LowerRLimitForTabulation = 1e-4;

        struct ZblParameters {
            std::array<Real, 6> coeffs;
            Real z1_z2;
            Real a_inverse;
        };

        Real cutoff;
        Real cutoff_transition_width;

        std::map<SpeciesSet<2, Symmetric>, ZblParameters> zbl_parameters;
        PerriotPolynomialCutoff cutoff_function;

        void loadDataset(std::istream& dataset, const std::set<SpeciesSet<2, Symmetric>>* species_filter = nullptr);

        // Builds the parameters for a pair, deducing z1_z2 and a_inverse from the species' atomic numbers.
        static ZblParameters makeParameters(const SpeciesSet<2, Symmetric>& pair, const std::array<Real, 6>& coeffs);

        std::array<Real, 2> energyAndDerivative(const ZblParameters& params, Real r) const;
    };
}

#endif