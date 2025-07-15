#include "core/potentials/ZblPotential.hpp"
#include "core/cutoff/PerriotPolynomialCutoff.hpp"

#include <fstream>

using namespace std;

namespace jgap {

    void ZblPotential::parseDmolFitCoefficients() {

        ifstream fIn(_dmolFile);
        nlohmann::json dmolFitCoefficients;
        fIn >> dmolFitCoefficients;

        _dmolFitCoefficients = {};
        for (auto& [key, val]: dmolFitCoefficients.items()) {
            auto Zs = split(key, ',');
            _dmolFitCoefficients[SpeciesPair{Z_inverse[stoi(Zs[0])], Z_inverse[stoi(Zs[1])]}] = {
                val[0], val[1], val[2], val[3], val[4], val[5]
            };
        }
    }

    ZblPotential::ZblPotential(const nlohmann::json& zblParams) {

        if (zblParams.contains("cutoff")) {
            _cutoff = zblParams["cutoff"]["cutoff"];
            _cutoffFunction = ParserRegistry<CutoffFunction>::get(zblParams["cutoff"]);
        } else {
            _cutoff = DEFAULT_ZBL_CUTOFF;
            _cutoffFunction = make_shared<PerriotPolynomialCutoff>(DEFAULT_ZBL_RMIN, _cutoff);
        }

        // TODO: not sure if default value is a good idea here

        _dmolFile = zblParams.value("dmol_fit_coefficients_file", DEFAULT_DMOL_FILE_PATH);
        parseDmolFitCoefficients();
    }

    nlohmann::json ZblPotential::serialize() {

        auto cutoffData = _cutoffFunction->serialize();
        cutoffData["type"] = _cutoffFunction->getType();

        if (_dmolFile != DEFAULT_DMOL_FILE_PATH) {
            return {
                {"cutoff", cutoffData},
                {"dmol_fit_coefficients_file", _dmolFile}
            };
        }

        return {
            {"cutoff", cutoffData}
        };
    }

    TabulationData ZblPotential::tabulate(const TabulationParams &params) {
        TabulationData result{};

        for (size_t i = 0; params.species.size(); i++) {
            for (size_t j = i; params.species.size(); j++) {
                auto speciesPair = SpeciesPair{params.species[i], params.species[j]};
                result.pairEnergies[speciesPair] = {};

                for (const double& r: params.grid2b) {
                    if (abs(r) < 1e-4) {
                        result.pairEnergies[speciesPair].push_back(0); // I don't like .eam.fs format :(
                    } else {
                        result.pairEnergies[speciesPair].push_back(
                            zblWithCutoff_eV(speciesPair, r)
                        );
                    }
                }
            }
        }

        return result;
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
