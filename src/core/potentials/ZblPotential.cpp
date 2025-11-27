#include "core/potentials/ZblPotential.hpp"
#include "core/cutoff/PerriotPolynomialCutoff.hpp"

#include <fstream>
#include <filesystem>

#if defined(_WIN32)
    #include <windows.h>
#elif defined(__APPLE__)
    #include <mach-o/dyld.h>
#else
    #include <unistd.h>
#endif


namespace jgap {

    ZblPotential::ZblPotential(const nlohmann::json& zblParams) {
        string test = zblParams.dump();
        JGAP_LOG_DEBUG("Parsing ZblPotential {}", zblParams.dump());
        if (zblParams.contains("cutoff")) {
            _cutoff = zblParams["cutoff"]["cutoff"];
            _cutoffFunction = ParserRegistry<CutoffFunction>::get(zblParams["cutoff"]);
        } else {
            _cutoff = DEFAULT_ZBL_CUTOFF;
            _cutoffFunction = make_shared<PerriotPolynomialCutoff>(DEFAULT_ZBL_R_MIN, DEFAULT_ZBL_CUTOFF);
        }

        _coeffFileName = zblParams.value("coefficients_file", "dmol-fit.json");

        ifstream fIn(_coeffFileName);
        if (!fIn.is_open()) {
            JGAP_LOG_WARN("Could not open coefficients_file: '" + _coeffFileName +
                                        "'. Trying to find it in resources.");
            fIn = ifstream(getResourcesCoeffFilePath(_coeffFileName));
            if (!fIn.is_open()) {
                JGAP_LOG_ERROR("Could not find coefficients_file in resources.", true);
            }
        }
        nlohmann::json dmolFitCoefficients;
        fIn >> dmolFitCoefficients;

        _dmolFitCoefficients = {};
        for (auto& [key, val]: dmolFitCoefficients.items()) {
            auto Zs = split(key, ',');
            _dmolFitCoefficients[SpeciesPair{CHEM_SYMBOLS[stoi(Zs[0])], CHEM_SYMBOLS[stoi(Zs[1])]}] = {
                val[0], val[1], val[2], val[3], val[4], val[5]
            };
        }

        _relevantSpecies = {};
        if (zblParams.contains("relevant_species")) {
            for (string species: zblParams["relevant_species"]) {
                _relevantSpecies.insert(species);
            }
        }
    }

    nlohmann::json ZblPotential::serialize() {

        auto cutoffData = _cutoffFunction->serialize();
        cutoffData["type"] = _cutoffFunction->getType();

        return {
            {"cutoff", cutoffData},
            {"coefficients_file", _coeffFileName},
            {"relevant_species", _relevantSpecies}
        };
    }

    void ZblPotential::tabulate(TabulationData &table) {

        for (Species species: table.allSpecies()) _relevantSpecies.insert(species);

        auto relevantSpeciesVector = vector(_relevantSpecies.begin(), _relevantSpecies.end());
        if (_relevantSpecies.empty()) {
            // TODO: tabulate all?
            JGAP_LOG_AND_THROW("Could not detect any relevant species to tabulate ZBL potential");
        }

        for (size_t i = 0; i < relevantSpeciesVector.size(); i++) {
            for (size_t j = i; j < relevantSpeciesVector.size(); j++) {
                auto speciesPair = SpeciesPair{relevantSpeciesVector[i], relevantSpeciesVector[j]};

                for (const auto& it: table.getOrMake2bGrid(speciesPair)) {
                    if (abs(it.pos) > _MINIMAL_TABULATED_R) {
                        it.value += zblWithCutoff_eV(speciesPair, it.pos);
                    }
                }
            }
        }
    }

    Predictions ZblPotential::predict(const AtomicStructure &structure) {

        double energy = 0;
        vector forces(structure.size(), Vector3{0, 0, 0});
        array<Vector3, 3> virials{};

        for (size_t i = 0; i < structure.size(); i++) {

            auto atom1 = structure[i];
            _relevantSpecies.insert(atom1.species());

            for (const NeighbourData& neighbour: atom1.neighbours()) {
                if (neighbour.index < i || neighbour.distance > _cutoff) continue;

                auto atom2 = structure[neighbour.index];

                double dE = zblWithCutoff_eV({atom1.species(), atom2.species()}, neighbour.distance);

                double dzbl_dr = zblWithCutoffDerivative_eV_per_Ang(
                    {atom1.species(), atom2.species()}, neighbour.distance
                    );
                Vector3 r21 = atom1.position() - (atom2.position() + neighbour.offset);
                Vector3 f21 = r21.normalize() * -dzbl_dr;

                if (neighbour.index == i) {
                    dE /= 2.0;
                    f21 /= 2.0;
                } else {
                    forces[i] = forces[i] + f21;
                    forces[neighbour.index] = forces[neighbour.index] - f21;
                }

                energy += dE;
                virials[0] += f21 * r21.x;
                virials[1] += f21 * r21.y;
                virials[2] += f21 * r21.z;
            }
        }

        return Predictions{
            energy,
            forces,
            virials
        };
    }

    string ZblPotential::getResourcesCoeffFilePath(const string& fileName) {
        filesystem::path exePath;

#if defined(_WIN32)
        char buffer[MAX_PATH];
        DWORD size = GetModuleFileNameA(NULL, buffer, MAX_PATH);
        if (size == 0 || size == MAX_PATH) {
            throw runtime_error("Failed to get executable path");
        }
        exePath = filesystem::path(buffer);

#elif defined(__APPLE__)
        char buffer[1024];
        uint32_t size = sizeof(buffer);
        if (_NSGetExecutablePath(buffer, &size) != 0) {
            throw runtime_error("Buffer too small for executable path");
        }
        exePath = filesystem::path(buffer);

#else // Linux and other Unix-like
        char buffer[1024];
        ssize_t count = readlink("/proc/self/exe", buffer, sizeof(buffer) - 1);
        if (count == -1) {
            throw std::runtime_error("Failed to read /proc/self/exe");
        }
        buffer[count] = '\0';
        exePath = filesystem::path(buffer);
#endif

        filesystem::path exeDir = exePath.parent_path();
        return (exeDir / "resources" / "dmol-screening-fit" / fileName).string();
    }

    double ZblPotential::zbl_eV(const SpeciesPair& speciesPair, double r) {
        auto coeffs = _dmolFitCoefficients[speciesPair];

        double Z1 = ATOMIC_NUMBERS[speciesPair.first()], Z2 = ATOMIC_NUMBERS[speciesPair.second()];
        double a = 0.46848 / (pow(Z1, 0.23) + pow(Z2, 0.23));

        double x = r / a;
        double phi = coeffs[0] * exp(-coeffs[1] * x)
                   + coeffs[2] * exp(-coeffs[3] * x)
                   + coeffs[4] * exp(-coeffs[5] * x);

        double r_meters = r * 1e-10;
        return Z1 * Z2 * _ELECTRON_CHARGE/*eV => no ^2*/ * phi / (4.0 * M_PI * _EPSILON * r_meters);
    }

    double ZblPotential::zblWithCutoff_eV(const SpeciesPair &speciesPair, const double r) {
        return _cutoffFunction->evaluate(r) * zbl_eV(speciesPair, r);
    }

    double ZblPotential::zblWithCutoffDerivative_eV_per_Ang(const SpeciesPair& speciesPair, const double r) {
        auto coeffs = _dmolFitCoefficients[speciesPair];

        double Z1 = ATOMIC_NUMBERS[speciesPair.first()], Z2 = ATOMIC_NUMBERS[speciesPair.second()];
        double a = 0.46848 / (pow(Z1, 0.23) + pow(Z2, 0.23));
        double r_meters = r * 1e-10;

        double x = r / a;
        double dx_dr = 1.0 / a;
        double phi = coeffs[0] * exp(-coeffs[1] * x)
                    + coeffs[2] * exp(-coeffs[3] * x)
                    + coeffs[4] * exp(-coeffs[5] * x);
        double dphi_dx = - coeffs[0] * coeffs[1] * exp(-coeffs[1] * x)
                         - coeffs[2] * coeffs[3] * exp(-coeffs[3] * x)
                         - coeffs[4] * coeffs[5] * exp(-coeffs[5] * x);
        double dphi_dr_meters = dx_dr * dphi_dx * 1e10;

        double dzbl_dr = - Z1 * Z2 * _ELECTRON_CHARGE/*eV => no ^2*/ * phi / (4.0 * M_PI * _EPSILON * r_meters * r_meters)
                         + Z1 * Z2 * _ELECTRON_CHARGE/*eV => no ^2*/ * dphi_dr_meters / (4.0 * M_PI * _EPSILON * r_meters);

        double dE_dr = _cutoffFunction->evaluate(r) * dzbl_dr * 1e-10
                       + _cutoffFunction->differentiate(r)/*already Angstrom*/ * zbl_eV(speciesPair, r);
        return dE_dr;
    }
}
