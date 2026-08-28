#include "ScreenedCoulombPotential.hpp"
#include "jgap/core/cutoff/PerriotPolynomialCutoff.hpp"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <sstream>

#include "ScreenedCoulombResources.h"
#include "jgap/core/atomic/iteration/Cluster2Expansion.hpp"
#include "jgap/core/atomic/neighbours/NeighbourLists.hpp"

namespace jgap {

    namespace {
        // Returns the contents of a built-in screening dataset. With #embed it is baked into the binary;
        // without it, the data is read at runtime from resources/dmol-screening-fit/<dataset>.dat, with a
        // clear error (and how to fix it) when that file can't be found.
        std::string getEmbeddedDataset(EmbeddedScreenedCoulombCoeffDataset type) {
#if defined(__has_embed)
            switch (type) {
                case EmbeddedScreenedCoulombCoeffDataset::DMOL:
                    return std::string(reinterpret_cast<const char*>(dmol_dat), dmol_dat_len);
                case EmbeddedScreenedCoulombCoeffDataset::MP2:
                    return std::string(reinterpret_cast<const char*>(mp2_dat), mp2_dat_len);
            }
            return "";
#else
            std::string dataset_name;
            switch (type) {
                case EmbeddedScreenedCoulombCoeffDataset::DMOL:
                    dataset_name = "dmol.dat";
                    break;
                case EmbeddedScreenedCoulombCoeffDataset::MP2:
                    dataset_name = "mp2.dat";
                    break;
            }
            const std::string path = "resources/dmol-screening-fit/" + dataset_name;

            std::ifstream file(path, std::ios::binary);
            if (!file.is_open()) {
                JGAP_LOG_AND_THROW(
                    "This build was compiled without #embed (needs GCC 15+ or Clang 19+), so the built-in "
                    "screening dataset is loaded from '{}' at runtime, but that file was not found. "
                    "Fix it by one of:\n"
                    "  - run from a directory that contains the project's 'resources/' folder (e.g. copy "
                    "resources/ next to the executable),\n"
                    "  - rebuild jgap with a #embed-capable compiler (GCC 15+ / Clang 19+), or\n"
                    "  - pass an explicit screening dataset file (ScreenedCoulombPotential(filename, training_data), "
                    "or the "
                    "fit's screened_coulomb_dataset_file option).",
                    path
                );
            }
            std::ostringstream contents;
            contents << file.rdbuf();
            return contents.str();
#endif
        }

        std::set<Species2Sorted> deducePairs(const std::vector<Atoms>& training_data) {
            std::set<Species2Sorted> species_pairs;
            for (const auto& atoms: training_data) {
                const auto& species_vec = atoms.getSpecies();
                std::set unique_species(species_vec.begin(), species_vec.end());
                for (auto it1 = unique_species.begin(); it1 != unique_species.end(); ++it1) {
                    for (auto it2 = it1; it2 != unique_species.end(); ++it2) {
                        species_pairs.insert(Species2Sorted{*it1, *it2});
                    }
                }
            }
            return species_pairs;
        }
    }

    void ScreenedCoulombPotential::loadDataset(std::istream& dataset, const std::set<Species2Sorted>* species_filter) {
        std::string line;
        while (std::getline(dataset, line)) {
            if (line.empty()) continue;
            std::istringstream iss(line);
            std::string elem1_str, elem2_str;
            std::array<Real, 6> coeffs = {};

            if (!(iss >> elem1_str >> elem2_str >> coeffs[0] >> coeffs[1] >> coeffs[2] >> coeffs[3] >> coeffs[4]
                  >> coeffs[5])) {
                continue;
            }

            Species s1(elem1_str);
            Species s2(elem2_str);
            Species2Sorted pair{s1, s2};

            if (species_filter && !species_filter->contains(pair)) {
                continue;
            }

            screened_coulomb_parameters[pair] = makeParameters(pair, coeffs);
        }
    }

    ScreenedCoulombPotential::ScreenedCoulombParameters ScreenedCoulombPotential::makeParameters(
        const Species2Sorted& pair, const std::array<Real, 6>& coeffs
    ) {
        const Species& s1 = pair.nodes[0];
        const Species& s2 = pair.nodes[1];

        if (!s1.atomicNumber().has_value()) {
            JGAP_LOG_AND_THROW("Atomic number undefined for {}", s1.symbol());
        }
        if (!s2.atomicNumber().has_value()) {
            JGAP_LOG_AND_THROW("Atomic number undefined for {}", s2.symbol());
        }

        Real z1 = s1.atomicNumber().value();
        Real z2 = s2.atomicNumber().value();

        return {coeffs, z1 * z2};
    }

    ScreenedCoulombPotential::ScreenedCoulombPotential(
        const std::set<Species2Sorted>& species,
        EmbeddedScreenedCoulombCoeffDataset embedded_dataset,
        Real cutoff,
        Real cutoff_transition_width
    ) :
        cutoff(cutoff),
        cutoff_transition_width(cutoff_transition_width),
        cutoff_function(cutoff, cutoff_transition_width) {
        std::istringstream dataset_stream(getEmbeddedDataset(embedded_dataset));
        loadDataset(dataset_stream, &species);
    }

    ScreenedCoulombPotential::ScreenedCoulombPotential(
        const std::vector<Atoms>& training_data,
        EmbeddedScreenedCoulombCoeffDataset embedded_dataset,
        Real cutoff,
        Real cutoff_transition_width
    ) :
        cutoff(cutoff),
        cutoff_transition_width(cutoff_transition_width),
        cutoff_function(cutoff, cutoff_transition_width) {
        std::set<Species2Sorted> pairs = deducePairs(training_data);
        std::istringstream dataset_stream(getEmbeddedDataset(embedded_dataset));
        loadDataset(dataset_stream, &pairs);
    }

    ScreenedCoulombPotential::ScreenedCoulombPotential(
        const std::string& dataset_filename,
        const std::vector<Atoms>& training_data,
        Real cutoff,
        Real cutoff_transition_width
    ) :
        cutoff(cutoff),
        cutoff_transition_width(cutoff_transition_width),
        cutoff_function(cutoff, cutoff_transition_width) {
        std::ifstream dataset_file(dataset_filename);
        if (!dataset_file.is_open()) {
            JGAP_LOG_AND_THROW("Could not open screening dataset file '{}'", dataset_filename);
        }
        std::set<Species2Sorted> pairs = deducePairs(training_data);
        loadDataset(dataset_file, &pairs);
    }

    ScreenedCoulombPotential::ScreenedCoulombPotential(
        std::istream& custom_dataset, const std::set<Species2Sorted>& species, Real cutoff, Real cutoff_transition_width
    ) :
        cutoff(cutoff),
        cutoff_transition_width(cutoff_transition_width),
        cutoff_function(cutoff, cutoff_transition_width) {
        loadDataset(custom_dataset, &species);
    }

    ScreenedCoulombPotential::ScreenedCoulombPotential(
        std::istream& custom_dataset, const std::vector<Atoms>& training_data, Real cutoff, Real cutoff_transition_width
    ) :
        cutoff(cutoff),
        cutoff_transition_width(cutoff_transition_width),
        cutoff_function(cutoff, cutoff_transition_width) {
        std::set<Species2Sorted> pairs = deducePairs(training_data);
        loadDataset(custom_dataset, &pairs);
    }

    ScreenedCoulombPotential::ScreenedCoulombPotential(
        const std::map<Species2Sorted, std::array<Real, 6>>& coefficients, Real cutoff, Real cutoff_transition_width
    ) :
        cutoff(cutoff),
        cutoff_transition_width(cutoff_transition_width),
        cutoff_function(cutoff, cutoff_transition_width) {
        for (const auto& [pair, coeffs]: coefficients) {
            screened_coulomb_parameters[pair] = makeParameters(pair, coeffs);
        }
    }

    std::map<Species2Sorted, std::array<Real, 6>> ScreenedCoulombPotential::getCoefficients() const {
        std::map<Species2Sorted, std::array<Real, 6>> coefficients;
        for (const auto& [pair, params]: screened_coulomb_parameters) {
            coefficients[pair] = params.coeffs;
        }
        return coefficients;
    }

    void ScreenedCoulombPotential::fillTables(TabulationData& tables) const {
        for (const auto& [species_pair, params]: screened_coulomb_parameters) {
            auto& table = tables.two_body_grids.getValueGrid(species_pair);
            for (const auto cell: table) {
                if (cell.pos[0] < LowerRLimitForTabulation) continue;
                cell.value += energyAndDerivative(params, cell.pos[0])[0];
            }
        }
    }

    AtomicQuantity ScreenedCoulombPotential::calculateEnergy(const Atoms& atoms) const {
        AtomicQuantity result(atoms.nAtoms());

        NeighbourLists nl(atoms, cutoff);

        for (const auto& [species_pair, params]: screened_coulomb_parameters) {
            Cluster2Expansion expansion(species_pair);
            expansion.forEach(nl, [&](const Cluster2& cluster) {
                auto [E_pair, dE_dr_pair] = energyAndDerivative(params, cluster.separation01.magnitude);

                Real E_cluster = 0.5 * E_pair;
                result.value += E_cluster;

                Real dE_dr = 0.5 * dE_dr_pair;
                utils::accumulatePairDistanceDerivatives(
                    result.forces[cluster.idx0],
                    result.forces[cluster.idx1],
                    result.virials,
                    dE_dr,
                    cluster.separation01
                );
            });
        }

        return result;
    }

    std::array<Real, 2> ScreenedCoulombPotential::energyAndDerivative(
        const Species2Sorted& species_pair, Real r
    ) const {
        auto it = screened_coulomb_parameters.find(species_pair);
        if (it == screened_coulomb_parameters.end()) {
            return {0.0, 0.0};
        }
        return energyAndDerivative(it->second, r);
    }

    std::array<Real, 2> ScreenedCoulombPotential::energyAndDerivative(
        const ScreenedCoulombParameters& params, Real r
    ) const {
        Real term1 = exp(-params.coeffs[1] * r);
        Real term2 = exp(-params.coeffs[3] * r);
        Real term3 = exp(-params.coeffs[5] * r);

        Real phi = params.coeffs[0] * term1 + params.coeffs[2] * term2 + params.coeffs[4] * term3;

        Real dphi_dr = -params.coeffs[0] * params.coeffs[1] * term1 - params.coeffs[2] * params.coeffs[3] * term2
                       - params.coeffs[4] * params.coeffs[5] * term3;

        Real prefactor = params.z1_z2 * CoulombConstant_eV_Ang / r;

        Real sc_e = prefactor * phi;
        Real dsc_dr = -(sc_e / r) + prefactor * dphi_dr;

        auto [cut, dcut] = cutoff_function.evaluateAndDifferentiate(r);

        Real E = cut * sc_e;
        Real dE_dr = cut * dsc_dr + dcut * sc_e;

        return {E, dE_dr};
    }
}
