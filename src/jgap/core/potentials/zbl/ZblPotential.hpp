#ifndef ZBLPOTENTIAL_HPP
#define ZBLPOTENTIAL_HPP

#include <array>
#include <istream>
#include <map>
#include <set>
#include <string>

#include "jgap/core/atomic/Atoms.hpp"
#include "jgap/core/atomic/species/composition/Species2Sorted.hpp"
#include "jgap/core/cutoff/CutoffFunction.hpp"
#include "jgap/core/cutoff/PerriotPolynomialCutoff.hpp"
#include "jgap/core/potentials/Potential.hpp"
#include "jgap/io/log/CurrentLogger.hpp"

namespace jgap {

    // The built-in screening datasets. When the compiler supports #embed they are baked into the binary;
    // otherwise they are read at runtime from resources/dmol-screening-fit/<dataset>.dat (see the .cpp).
    enum class EmbeddedZBLCoeffDataset { DMOL, NLH, MP2 };

    class ZblPotential : public Potential {
    public:
        static constexpr Real DefaultZblCutoff = 2.2;
        static constexpr Real DefaultZblCutoffTransitionWidth = 1.2;
        static constexpr Real Epsilon0_F_per_m = 8.854187817e-12;
        static constexpr Real ElectronCharge_C = 1.60217657e-19;
        static constexpr Real CoulombConstant_eV_Ang = ElectronCharge_C / (4.0 * M_PI * Epsilon0_F_per_m * 1e-10);

        // Built-in dataset (embedded via #embed, or read from resources/ at runtime when #embed is
        // unavailable; see the .cpp), restricted to the given pairs / to the elements in the training data.
        ZblPotential(const std::set<Species2Sorted> &species,
                     EmbeddedZBLCoeffDataset embedded_dataset = EmbeddedZBLCoeffDataset::DMOL,
                     Real cutoff = DefaultZblCutoff, Real cutoff_transition_width = DefaultZblCutoffTransitionWidth);

        ZblPotential(const std::vector<Atoms> &training_data,
                     EmbeddedZBLCoeffDataset embedded_dataset = EmbeddedZBLCoeffDataset::DMOL,
                     Real cutoff = DefaultZblCutoff, Real cutoff_transition_width = DefaultZblCutoffTransitionWidth);

        ZblPotential(std::istream &custom_dataset, const std::set<Species2Sorted> &species,
                     Real cutoff = DefaultZblCutoff, Real cutoff_transition_width = DefaultZblCutoffTransitionWidth);

        ZblPotential(std::istream &custom_dataset, const std::vector<Atoms> &training_data,
                     Real cutoff = DefaultZblCutoff, Real cutoff_transition_width = DefaultZblCutoffTransitionWidth);

        // Reads coefficients from a dataset file, keeping only the element pairs present in the training
        // data.
        ZblPotential(const std::string &dataset_filename, const std::vector<Atoms> &training_data,
                     Real cutoff = DefaultZblCutoff, Real cutoff_transition_width = DefaultZblCutoffTransitionWidth);

        // z1_z2 and a_inverse are deduced from the species in each pair.
        ZblPotential(const std::map<Species2Sorted, std::array<Real, 6>> &coefficients, Real cutoff = DefaultZblCutoff,
                     Real cutoff_transition_width = DefaultZblCutoffTransitionWidth);

        std::map<Species2Sorted, std::array<Real, 6>> getCoefficients() const;

        Real getCutoff() const { return cutoff; }
        Real getCutoffTransitionWidth() const { return cutoff_transition_width; }

        AtomicQuantity calculateEnergy(const Atoms &atoms) const override;

        std::array<Real, 2> energyAndDerivative(const Species2Sorted &species_pair, Real r) const;

        Cutoffs getCutoffs() const override { return {{2u, cutoff}}; }

        void fillTables(TabulationData &tables) const override;

        ZblPotential *clone() const override { return new ZblPotential(*this); }

    private:
        static constexpr Real LowerRLimitForTabulation = 1e-4;

        struct ZblParameters {
            std::array<Real, 6> coeffs;
            Real z1_z2;
            Real a_inverse;
        };

        Real cutoff;
        Real cutoff_transition_width;

        std::map<Species2Sorted, ZblParameters> zbl_parameters;
        PerriotPolynomialCutoff cutoff_function;

        void loadDataset(std::istream &dataset, const std::set<Species2Sorted> *species_filter = nullptr);

        // Builds the parameters for a pair, deducing z1_z2 and a_inverse from the species' atomic numbers.
        static ZblParameters makeParameters(const Species2Sorted &pair, const std::array<Real, 6> &coeffs);

        std::array<Real, 2> energyAndDerivative(const ZblParameters &params, Real r) const;
    };
}

#endif
