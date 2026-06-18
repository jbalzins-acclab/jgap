#include <gtest/gtest.h>
#include <core/potentials/zbl/ZblPotential.hpp>
#include "core/atomic/species/Species.hpp"
#include "core/atomic/Atoms.hpp"

#include <sstream>
#include <iostream>
#include <vector>
#include <array>
#include <cmath>

using namespace std;
using namespace jgap;

TEST(ZblPotentialTest, FeNiInteractionFromEmbedded) {
    Species fe("Fe"), ni("Ni");
    std::set species = {SpeciesSet<2, Symmetric>{fe, ni}};

    ZblPotential zbl(species, EmbeddedZBLCoeffDataset::DMOL, 10.0);

    // Test at a specific distance
    Real r = 0.5; // Angstroms
    auto [energy, derivative] = zbl.energyAndDerivative(SpeciesSet<2, Symmetric>{fe, ni}, r);

    // The coefficients for Fe-Ni from dmol.dat are:
    // 0.32818750683529224 1.2763590338471473 0.5319935400634741 0.4920242580140316 0.13562223656825356 2.9359601944688674

    Real z1 = fe.atomicNumber().value();
    Real z2 = ni.atomicNumber().value();
    std::array c = {
        0.32818750683529224, 1.2763590338471473,
        0.5319935400634741,  0.4920242580140316,
        0.13562223656825356, 2.9359601944688674
    };

    Real a_inverse = (std::pow(z1, 0.23) + std::pow(z2, 0.23)) / 0.46848;
    Real x = r * a_inverse;

    Real term1 = exp(-c[1] * x);
    Real term2 = exp(-c[3] * x);
    Real term3 = exp(-c[5] * x);

    Real phi = c[0] * term1 + c[2] * term2 + c[4] * term3;
    Real dphi_dx = -c[0]*c[1]*term1 - c[2]*c[3]*term2 - c[4]*c[5]*term3;
    Real dphi_dr = a_inverse * dphi_dx;

    Real prefactor = z1 * z2 * ZblPotential::CoulombConstant_eV_Ang / r;
    Real expected_e = prefactor * phi;
    Real expected_dE_dr = -(prefactor * phi / r) + prefactor * dphi_dr;

    ASSERT_NEAR(1206.6878058678733, expected_e, 1e-6); // dmol-fit.py prediction
    ASSERT_NEAR(energy, expected_e, 1e-6);
    ASSERT_NEAR(derivative, expected_dE_dr, 1e-6);
}

TEST(ZblPotentialTest, MockInteraction) {
    std::stringstream mock_stream;
    mock_stream << "H He 0.1 0.2 0.3 0.4 0.5 0.6" << std::endl;

    Species h("H"), he("He");
    std::set species = {SpeciesSet<2, Symmetric>{h, he}};

    ZblPotential zbl(mock_stream, species, 10.0);

    Real r = 1.0;
    auto [energy, derivative] = zbl.energyAndDerivative(SpeciesSet<2, Symmetric>{h, he}, r);

    // --- Manual Calculation ---
    // Z(H) = 1, Z(He) = 2
    // a_inverse = (1^0.23 + 2^0.23) / 0.46848 = (1 + 1.173) / 0.46848 = 4.638
    // x = 1.0 * 4.638 = 4.638
    // phi = 0.1*exp(-0.2*x) + 0.3*exp(-0.4*x) + 0.5*exp(-0.6*x)
    // phi = 0.1*exp(-0.9276) + 0.3*exp(-1.8552) + 0.5*exp(-2.7828)
    // phi = 0.1*0.3955 + 0.3*0.1564 + 0.5*0.0619 = 0.03955 + 0.04692 + 0.03095 = 0.11742
    // zbl_e = (1*2*14.3996/1.0) * 0.11742 = 3.3816

    // dphi_dx = -0.1*0.2*exp(-0.9276) - 0.3*0.4*exp(-1.8552) - 0.5*0.6*exp(-2.7828)
    // dphi_dx = -0.02*0.3955 - 0.12*0.1564 - 0.3*0.0619 = -0.00791 - 0.018768 - 0.01857 = -0.045248
    // dphi_dr = 4.638 * -0.045248 = -0.2100
    // dzbl_dr = -(3.3816/1.0) + (28.7992/1.0) * -0.2100 = -3.3816 - 6.0478 = -9.4294

    ASSERT_NEAR(energy, 3.3816, 1e-3);
    ASSERT_NEAR(derivative, -9.4294, 1e-2);
}

TEST(ZblPotentialTest, Mock2AtomSystem) {
    std::stringstream mock_stream;
    mock_stream << "H He 0.1 0.2 0.3 0.4 0.5 0.6\n";

    Species h("H"), he("He");
    std::set species = {SpeciesSet<2, Symmetric>{h, he}};

    ZblPotential zbl(mock_stream, species);

    Atoms atoms({ {0, 0, 0}, {1.0, 0, 0} }, {h, he});
    AtomicQuantity quantity = zbl.calculateEnergy(atoms);

    auto [expected_energy, expected_deriv] = zbl.energyAndDerivative(SpeciesSet<2, Symmetric>{h, he}, 1.0);

    // Check energy
    EXPECT_NEAR(quantity.value, expected_energy, 1e-6);

    // Forces: dE_dr is negative (repulsive).
    // H is at origin, He is at +x.
    // H gets pushed to -x, He gets pushed to +x.
    // The magnitude of the force is -dE_dr.
    EXPECT_NEAR(quantity.forces[0].x, expected_deriv, 1e-6);
    EXPECT_NEAR(quantity.forces[1].x, -expected_deriv, 1e-6);

    // Virials check
    EXPECT_NEAR(quantity.virials.xx, -1.0 * expected_deriv, 1e-6);
}

TEST(ZblPotentialTest, Mock3AtomSystem) {
    std::stringstream mock_stream;
    mock_stream << "H He 0.1 0.2 0.3 0.4 0.5 0.6\n";
    mock_stream << "H Li 0.2 0.3 0.4 0.5 0.6 0.7\n";
    mock_stream << "He Li 0.3 0.4 0.5 0.6 0.7 0.8\n";

    Species h("H"), he("He"), li("Li");
    std::set species = {SpeciesSet<2, Symmetric>{h, he}, SpeciesSet<2, Symmetric>{h, li}, SpeciesSet<2, Symmetric>{he, li}};

    ZblPotential zbl(mock_stream, species);

    Atoms atoms(
        { {0, 0, 0}, {1.0, 0, 0}, {0, 1.0, 0} },
        {h, he, li}
    );

    AtomicQuantity quantity = zbl.calculateEnergy(atoms);

    auto [e_h_he, dE_h_he] = zbl.energyAndDerivative(SpeciesSet<2, Symmetric>{h, he}, 1.0);
    auto [e_h_li, dE_h_li] = zbl.energyAndDerivative(SpeciesSet<2, Symmetric>{h, li}, 1.0);
    Real r_he_li = std::sqrt(2.0);
    auto [e_he_li, dE_he_li] = zbl.energyAndDerivative(SpeciesSet<2, Symmetric>{he, li}, r_he_li);

    // Check total energy
    Real expected_total_energy = e_h_he + e_h_li + e_he_li;
    EXPECT_NEAR(quantity.value, expected_total_energy, 1e-6);

    // Check forces
    // H (0) is repelled by He (1) along -x, by Li (2) along -y
    EXPECT_NEAR(quantity.forces[0].x, dE_h_he, 1e-6);
    EXPECT_NEAR(quantity.forces[0].y, dE_h_li, 1e-6);
    EXPECT_NEAR(quantity.forces[0].z, 0.0, 1e-6);

    // He (1) is repelled by H (0) along +x
    // He (1) interacts with Li (2). Vector from Li to He is (1, -1, 0).
    // The force on He due to Li is -dE_he_li * (1/sqrt(2), -1/sqrt(2), 0)
    Real f_he_li_mag = -dE_he_li;
    EXPECT_NEAR(quantity.forces[1].x, -dE_h_he + f_he_li_mag * (1.0 / r_he_li), 1e-6);
    EXPECT_NEAR(quantity.forces[1].y, f_he_li_mag * (-1.0 / r_he_li), 1e-6);

    // Li (2) is repelled by H (0) along +y
    // Li (2) interacts with He (1). The force on Li due to He is -dE_he_li * (-1/sqrt(2), 1/sqrt(2), 0)
    EXPECT_NEAR(quantity.forces[2].x, f_he_li_mag * (-1.0 / r_he_li), 1e-6);
    EXPECT_NEAR(quantity.forces[2].y, -dE_h_li + f_he_li_mag * (1.0 / r_he_li), 1e-6);

    // Check sum of forces is zero (Newton's third law)
    Vector3 force_sum = quantity.forces[0] + quantity.forces[1] + quantity.forces[2];
    EXPECT_NEAR(force_sum.x, 0.0, 1e-6);
    EXPECT_NEAR(force_sum.y, 0.0, 1e-6);
    EXPECT_NEAR(force_sum.z, 0.0, 1e-6);
}