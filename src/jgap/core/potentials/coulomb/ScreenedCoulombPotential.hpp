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
        static constexpr double DefaultCutoff = 2.2;
        static constexpr double DefaultCutoffTransitionWidth = 1.2;
        static constexpr double Epsilon0_F_per_m = 8.854187817e-12;
        static constexpr double ElectronCharge_C = 1.60217657e-19;
        static constexpr double CoulombConstant_eV_Ang = ElectronCharge_C / (4.0 * M_PI * Epsilon0_F_per_m * 1e-10);

        // Built-in dataset (embedded via #embed, or read from resources/ at runtime when #embed is
        // unavailable; see the .cpp), restricted to the given pairs / to the elements in the training data.
        ScreenedCoulombPotential(
            const std::set<Species2Sorted>& species,
            EmbeddedScreenedCoulombCoeffDataset embedded_dataset = EmbeddedScreenedCoulombCoeffDataset::DMOL, double cutoff = DefaultCutoff,
            double cutoff_transition_width = DefaultCutoffTransitionWidth
        );

        ScreenedCoulombPotential(
            const std::vector<Atoms>& training_data,
            EmbeddedScreenedCoulombCoeffDataset embedded_dataset = EmbeddedScreenedCoulombCoeffDataset::DMOL, double cutoff = DefaultCutoff,
            double cutoff_transition_width = DefaultCutoffTransitionWidth
        );

        ScreenedCoulombPotential(
            std::istream& custom_dataset, const std::set<Species2Sorted>& species, double cutoff = DefaultCutoff,
            double cutoff_transition_width = DefaultCutoffTransitionWidth
        );

        ScreenedCoulombPotential(
            std::istream& custom_dataset, const std::vector<Atoms>& training_data, double cutoff = DefaultCutoff,
            double cutoff_transition_width = DefaultCutoffTransitionWidth
        );

        // Reads coefficients from a dataset file, keeping only the element pairs present in the training
        // data.
        ScreenedCoulombPotential(
            const std::string& dataset_filename, const std::vector<Atoms>& training_data,
            double cutoff = DefaultCutoff, double cutoff_transition_width = DefaultCutoffTransitionWidth
        );

        // z1_z2 and a_inverse are deduced from the species in each pair.
        ScreenedCoulombPotential(
            const std::map<Species2Sorted, std::array<double, 6>>& coefficients, double cutoff = DefaultCutoff,
            double cutoff_transition_width = DefaultCutoffTransitionWidth
        );

        std::map<Species2Sorted, std::array<double, 6>> getCoefficients() const;

        double getCutoff() const { return cutoff; }
        double getCutoffTransitionWidth() const { return cutoff_transition_width; }

        AtomicQuantity calculateEnergy(const Atoms& atoms) const override;

        std::array<double, 2> energyAndDerivative(const Species2Sorted& species_pair, double r) const;

        Cutoffs getCutoffs() const override { return {{2u, cutoff}}; }

        void fillTables(TabulationData& tables) const override;

        ScreenedCoulombPotential* clone() const override { return new ScreenedCoulombPotential(*this); }

    private:
        static constexpr double LowerRLimitForTabulation = 1e-4;

        struct ScreenedCoulombParameters {
            std::array<double, 6> coeffs;
            double z1_z2;
        };

        double cutoff;
        double cutoff_transition_width;

        std::map<Species2Sorted, ScreenedCoulombParameters> screened_coulomb_parameters;
        PerriotPolynomialCutoff cutoff_function;

        void loadDataset(std::istream& dataset, const std::set<Species2Sorted>* species_filter = nullptr);

        // Builds the parameters for a pair, deducing z1_z2 from the species' atomic numbers.
        static ScreenedCoulombParameters makeParameters(const Species2Sorted& pair, const std::array<double, 6>& coeffs);

        std::array<double, 2> energyAndDerivative(const ScreenedCoulombParameters& params, double r) const;
    };
}

#endif
