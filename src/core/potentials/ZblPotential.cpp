#include "core/potentials/ZblPotential.hpp"

#include "core/cutoff/PerriotPolynomialCutoff.hpp"

namespace jgap {

    ZblPotential::ZblPotential(nlohmann::json dmolFitCoefficients, double cutoffLow, double cutoffHigh) {

        StdoutLogger::initIfNotInitialized();

        _cutoff = cutoffHigh;
        _cutoffFunction = make_shared<PerriotPolynomialCutoff>(cutoffLow, cutoffHigh);

        _dmolFitCoefficients = {};
        for (auto& [key, val]: dmolFitCoefficients.items()) {
            auto Zs = split(key, ',');
            _dmolFitCoefficients[SpeciesPair{Z_inverse[stoi(Zs[0])], Z_inverse[stoi(Zs[1])]}] = {
                val[0], val[1], val[2], val[3], val[4], val[5]
            };
        }
    }

    PotentialPrediction ZblPotential::predict(const AtomicStructure &structure) {

        double energy = 0;
        vector forces(structure.atoms.size(), Vector3{0, 0, 0});
        // TODO: virials?
        for (size_t i = 0; i < structure.atoms.size(); i++) {

            auto atom1 = structure.atoms[i];

            for (const NeighbourData& neighbour: atom1.neighbours.value()) {
                if (neighbour.index < i || neighbour.distance > _cutoff) continue;

                auto atom2 = structure.atoms[neighbour.index];

                double dE = zblWithCutoff_eV({atom1.species, atom2.species}, neighbour.distance);

                if (neighbour.index == i) {
                    energy += dE / 2;
                    continue;
                }

                energy += dE;

                double dzbl_dr = zblWithCutoffDerivative_eV_per_Ang({atom1.species, atom2.species}, neighbour.distance);
                Vector3 f = (atom1.position - (atom2.position + neighbour.offset)).normalize() * dzbl_dr;
                forces[i] = forces[i] - f;
                forces[neighbour.index] = forces[neighbour.index] + f;
            }
        }

        return PotentialPrediction{
            energy,
            forces
        };
    }

    nlohmann::json ZblPotential::serialize() {
        return {"name", "zbl"};
    }

    double ZblPotential::zbl_eV(const SpeciesPair& speciesPair, double r) {
        auto coeffs = _dmolFitCoefficients[speciesPair];

        /*
        def screened_coulomb_ev(r, c1, c2, c3, c4, c5, c6):
        eps = 8.854187817e-12
        e = 1.60217657e-19

        a = 0.46848 / (Z1**0.23 + Z2**0.23)
        x = r / a
        phi = c1 * np.exp(-c2 * x) + c3 * np.exp(-c4 * x) + c5 * np.exp(-c6 * x)
        r = r * 1e-10  # Å --> m
        E = Z1 * Z2 * e * phi / (4.0 * np.pi * eps * r) # not e**2 => eV
        return E
         */

        double Z1 = Z_default[speciesPair.first()], Z2 = Z_default[speciesPair.second()];
        double a = 0.46848 / (pow(Z1, 0.23) + pow(Z2, 0.23));

        double x = r / a;
        double phi = coeffs[0] * exp(-coeffs[1] * x)
                    + coeffs[2] * exp(-coeffs[3] * x)
                    + coeffs[4] * exp(-coeffs[5] * x);

        double r_meters = r * 1e-10;
        return Z1 * Z2 * _electronCharge * phi / (4.0 * M_PI * _eps * r_meters);
    }

    double ZblPotential::zblWithCutoff_eV(const SpeciesPair &speciesPair, double r) {
        return _cutoffFunction->evaluate(r) * zbl_eV(speciesPair, r);
    }

    double ZblPotential::zblWithCutoffDerivative_eV_per_Ang(const SpeciesPair& speciesPair, double r) {
        auto coeffs = _dmolFitCoefficients[speciesPair];

        double Z1 = Z_default[speciesPair.first()], Z2 = Z_default[speciesPair.second()];
        double a = 0.46848 / (pow(Z1, 0.23) + pow(Z2, 0.23));
        double r_meters = r * 1e-10;

        double x = r / a;
        double dx_dr = 1 / a;
        double phi = coeffs[0] * exp(-coeffs[1] * x)
                    + coeffs[2] * exp(-coeffs[3] * x)
                    + coeffs[4] * exp(-coeffs[5] * x);
        double dphi_dx = - coeffs[0] * coeffs[1] * exp(-coeffs[1] * x)
                         - coeffs[2] * coeffs[3] * exp(-coeffs[3] * x)
                         - coeffs[4] * coeffs[5] * exp(-coeffs[5] * x);
        double dphi_drmeters = dx_dr * dphi_dx * 1e10;

        double dzbl_dr = - Z1 * Z2 * _electronCharge * phi / (4.0 * M_PI * _eps * r_meters * r_meters)
                         + Z1 * Z2 * _electronCharge * dphi_drmeters / (4.0 * M_PI * _eps * r_meters);

        double dE_dr = _cutoffFunction->evaluate(r) * dzbl_dr * 1e-10
                       + _cutoffFunction->differentiate(r)/*already Angstrom*/ * zbl_eV(speciesPair, r);
        return dE_dr;
    }
}
