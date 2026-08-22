#include <gtest/gtest.h>
#include <jgap/core/potentials/coulomb/ScreenedCoulombPotential.hpp>
#include "jgap/core/atomic/Atoms.hpp"
#include "jgap/core/atomic/species/Species.hpp"
#include "jgap/core/atomic/species/composition/Species2Sorted.hpp"

#include <array>
#include <cmath>

using namespace std;
using namespace jgap;

TEST(ScreenedCoulombPotentialTest, FeNiInteractionFromEmbedded) {
    Species fe("Fe"), ni("Ni");
    std::set species = {Species2Sorted(fe, ni)};

    ScreenedCoulombPotential sc(species, EmbeddedScreenedCoulombCoeffDataset::DMOL, 10.0);

    // Test at a specific distance
    Real r = 0.5; // Angstroms
    auto [energy, derivative] = sc.energyAndDerivative(Species2Sorted(fe, ni), r);

    Real z1 = fe.atomicNumber().value();
    Real z2 = ni.atomicNumber().value();

    // The coefficients for Fe-Ni obtained from an old fit-dmol.py which uses alternative calculation
    // (with constant a=(std::pow(z1, 0.23) + std::pow(z2, 0.23)) / 0.46848 to make screening argument dimensionless):
    std::array c = {0.32818750683529224, 1.2763590338471473,  0.5319935400634741,
                    0.4920242580140316,  0.13562223656825356, 2.9359601944688674};

    Real a_inverse = (std::pow(z1, 0.23) + std::pow(z2, 0.23)) / 0.46848;
    Real x = r * a_inverse;

    Real term1 = exp(-c[1] * x);
    Real term2 = exp(-c[3] * x);
    Real term3 = exp(-c[5] * x);

    Real phi = c[0] * term1 + c[2] * term2 + c[4] * term3;
    Real dphi_dx = -c[0] * c[1] * term1 - c[2] * c[3] * term2 - c[4] * c[5] * term3;
    Real dphi_dr = a_inverse * dphi_dx;

    Real prefactor = z1 * z2 * ScreenedCoulombPotential::CoulombConstant_eV_Ang / r;
    Real expected_e = prefactor * phi;
    Real expected_dE_dr = -(prefactor * phi / r) + prefactor * dphi_dr;

    ASSERT_NEAR(1206.6878058678733, expected_e, 2e-3);
    ASSERT_NEAR(energy, 1206.6878058678733, 2e-3);
    ASSERT_NEAR(-7968.86, expected_dE_dr, 1e-2);
    ASSERT_NEAR(derivative, -7968.86, 2e-2);
}

TEST(ScreenedCoulombPotentialTest, FeNi4AtomSquare) {
    Species fe("Fe"), ni("Ni");
    std::set species = {Species2Sorted(fe, ni)};

    ScreenedCoulombPotential sc(species, EmbeddedScreenedCoulombCoeffDataset::DMOL, 10.0);

    Real r = 0.5; // side length of square in Angstroms
    auto [energy_single, deriv_single] = sc.energyAndDerivative(Species2Sorted(fe, ni), r);

    // 4 atoms in a square of side 0.5 Å:
    // Atom 0 (Fe): (0.0, 0.0, 0.0)
    // Atom 1 (Ni): (0.5, 0.0, 0.0)
    // Atom 2 (Fe): (0.5, 0.5, 0.0)
    // Atom 3 (Ni): (0.0, 0.5, 0.0)
    // Same species at opposing corners (Fe at 0 & 2, Ni at 1 & 3). Only Fe-Ni interaction is enabled.
    Atoms atoms({{0.0, 0.0, 0.0}, {r, 0.0, 0.0}, {r, r, 0.0}, {0.0, r, 0.0}}, {fe, ni, fe, ni});

    AtomicQuantity quantity = sc.calculateEnergy(atoms);

    // Total energy is 4x the single-pair energy
    EXPECT_NEAR(quantity.value, 4.0 * energy_single, 1e-5);

    // Forces:
    // Atom 0 at (0, 0): repelled by Atom 1 at (r, 0) and Atom 3 at (0, r)
    EXPECT_NEAR(quantity.forces[0].x, deriv_single, 1e-4);
    EXPECT_NEAR(quantity.forces[0].y, deriv_single, 1e-4);
    EXPECT_NEAR(quantity.forces[0].z, 0.0, 1e-6);

    // Atom 1 at (r, 0): repelled by Atom 0 at (0, 0) and Atom 2 at (r, r)
    EXPECT_NEAR(quantity.forces[1].x, -deriv_single, 1e-4);
    EXPECT_NEAR(quantity.forces[1].y, deriv_single, 1e-4);

    // Atom 2 at (r, r): repelled by Atom 1 at (r, 0) and Atom 3 at (0, r)
    EXPECT_NEAR(quantity.forces[2].x, -deriv_single, 1e-4);
    EXPECT_NEAR(quantity.forces[2].y, -deriv_single, 1e-4);

    // Atom 3 at (0, r): repelled by Atom 0 at (0, 0) and Atom 2 at (r, r)
    EXPECT_NEAR(quantity.forces[3].x, deriv_single, 1e-4);
    EXPECT_NEAR(quantity.forces[3].y, -deriv_single, 1e-4);
}
