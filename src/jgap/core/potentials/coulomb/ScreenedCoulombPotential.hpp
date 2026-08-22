#ifndef SCREENEDCOULOMBPOTENTIAL_HPP
#define SCREENEDCOULOMBPOTENTIAL_HPP

#include <array>
#include <istream>
#include <map>
#include <set>
#include <string>

#include "jgap/core/io/log/CurrentLogger.hpp"
#include "jgap/core/atomic/Atoms.hpp"
#include "jgap/core/atomic/species/composition/Species2Sorted.hpp"
#include "jgap/core/cutoff/CutoffFunction.hpp"
#include "jgap/core/cutoff/PerriotPolynomialCutoff.hpp"
#include "jgap/core/potentials/Potential.hpp"

namespace jgap {

    // The built-in screening datasets. When the compiler supports #embed they are baked into the binary;
    // otherwise they are read at runtime from resources/dmol-screening-fit/<dataset>.dat (see the .cpp).
    enum class EmbeddedScreenedCoulombCoeffDataset { DMOL, MP2 };

    class ScreenedCoulombPotential : public Potential {
    public:
        static constexpr Real DefaultCutoff = 2.2;
        static constexpr Real DefaultCutoffTransitionWidth = 1.2;
        static constexpr Real Epsilon0_F_per_m = 8.854187817e-12;
        static constexpr Real ElectronCharge_C = 1.60217657e-19;
        static constexpr Real CoulombConstant_eV_Ang = ElectronCharge_C / (4.0 * M_PI * Epsilon0_F_per_m * 1e-10);

        // Built-in dataset (embedded via #embed, or read from resources/ at runtime when #embed is
        // unavailable; see the .cpp), restricted to the given pairs / to the elements in the training data.
        ScreenedCoulombPotential(
            const std::set<Species2Sorted>& species,
            EmbeddedScreenedCoulombCoeffDataset embedded_dataset = EmbeddedScreenedCoulombCoeffDataset::DMOL, Real cutoff = DefaultCutoff,
            Real cutoff_transition_width = DefaultCutoffTransitionWidth
        );

        ScreenedCoulombPotential(
            const std::vector<Atoms>& training_data,
            EmbeddedScreenedCoulombCoeffDataset embedded_dataset = EmbeddedScreenedCoulombCoeffDataset::DMOL, Real cutoff = DefaultCutoff,
            Real cutoff_transition_width = DefaultCutoffTransitionWidth
        );

        ScreenedCoulombPotential(
            std::istream& custom_dataset, const std::set<Species2Sorted>& species, Real cutoff = DefaultCutoff,
            Real cutoff_transition_width = DefaultCutoffTransitionWidth
        );

        ScreenedCoulombPotential(
            std::istream& custom_dataset, const std::vector<Atoms>& training_data, Real cutoff = DefaultCutoff,
            Real cutoff_transition_width = DefaultCutoffTransitionWidth
        );

        // Reads coefficients from a dataset file, keeping only the element pairs present in the training
        // data.
        ScreenedCoulombPotential(
            const std::string& dataset_filename, const std::vector<Atoms>& training_data,
            Real cutoff = DefaultCutoff, Real cutoff_transition_width = DefaultCutoffTransitionWidth
        );

        // z1_z2 and a_inverse are deduced from the species in each pair.
        ScreenedCoulombPotential(
            const std::map<Species2Sorted, std::array<Real, 6>>& coefficients, Real cutoff = DefaultCutoff,
            Real cutoff_transition_width = DefaultCutoffTransitionWidth
        );

        std::map<Species2Sorted, std::array<Real, 6>> getCoefficients() const;

        Real getCutoff() const { return cutoff; }
        Real getCutoffTransitionWidth() const { return cutoff_transition_width; }

        AtomicQuantity calculateEnergy(const Atoms& atoms) const override;

        std::array<Real, 2> energyAndDerivative(const Species2Sorted& species_pair, Real r) const;

        Cutoffs getCutoffs() const override { return {{2u, cutoff}}; }

        void fillTables(TabulationData& tables) const override;

        ScreenedCoulombPotential* clone() const override { return new ScreenedCoulombPotential(*this); }

    private:
        static constexpr Real LowerRLimitForTabulation = 1e-4;

        struct ScreenedCoulombParameters {
            std::array<Real, 6> coeffs;
            Real z1_z2;
        };

        Real cutoff;
        Real cutoff_transition_width;

        std::map<Species2Sorted, ScreenedCoulombParameters> screened_coulomb_parameters;
        PerriotPolynomialCutoff cutoff_function;

        void loadDataset(std::istream& dataset, const std::set<Species2Sorted>* species_filter = nullptr);

        // Builds the parameters for a pair, deducing z1_z2 from the species' atomic numbers.
        static ScreenedCoulombParameters makeParameters(const Species2Sorted& pair, const std::array<Real, 6>& coeffs);

        std::array<Real, 2> energyAndDerivative(const ScreenedCoulombParameters& params, Real r) const;
    };
}

#endif
