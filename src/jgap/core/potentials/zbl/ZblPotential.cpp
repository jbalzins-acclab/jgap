#include "ZblPotential.hpp"
#include "jgap/core/cutoff/PerriotPolynomialCutoff.hpp"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <sstream>

#include "ZblResources.h"
#include "jgap/core/atomic/iteration/Cluster2Expansion.hpp"
#include "jgap/core/atomic/neighbours/NeighbourLists.hpp"

namespace jgap {

    namespace {
        // Returns the contents of a built-in screening dataset. With #embed it is baked into the binary;
        // without it, the data is read at runtime from resources/dmol-screening-fit/<dataset>.dat, with a
        // clear error (and how to fix it) when that file can't be found.
        std::string getEmbeddedDataset(EmbeddedZBLCoeffDataset type) {
#if defined(__has_embed)
            switch (type) {
                case EmbeddedZBLCoeffDataset::NLH:
                    return std::string(reinterpret_cast<const char*>(nlh_dat), nlh_dat_len);
                case EmbeddedZBLCoeffDataset::DMOL:
                    return std::string(reinterpret_cast<const char*>(dmol_dat), dmol_dat_len);
                case EmbeddedZBLCoeffDataset::MP2:
                    return std::string(reinterpret_cast<const char*>(mp2_dat), mp2_dat_len);
            }
            return "";
#else
            std::string dataset_name;
            switch (type) {
                case EmbeddedZBLCoeffDataset::NLH:
                    dataset_name = "nlh.dat";
                    break;
                case EmbeddedZBLCoeffDataset::DMOL:
                    dataset_name = "dmol.dat";
                    break;
                case EmbeddedZBLCoeffDataset::MP2:
                    dataset_name = "mp2.dat";
                    break;
            }
            const std::string path = "resources/dmol-screening-fit/" + dataset_name;

            std::ifstream file(path, std::ios::binary);
            if (!file.is_open()) {
                JGAP_LOG_AND_THROW(
                    "This build was compiled without #embed (needs GCC 15+ or Clang 19+), so the built-in "
                    "ZBL screening dataset is loaded from '{}' at runtime, but that file was not found. "
                    "Fix it by one of:\n"
                    "  - run from a directory that contains the project's 'resources/' folder (e.g. copy "
                    "resources/ next to the executable),\n"
                    "  - rebuild jgap with a #embed-capable compiler (GCC 15+ / Clang 19+), or\n"
                    "  - pass an explicit ZBL dataset file (ZblPotential(filename, training_data), or the "
                    "fit's zbl_dataset_file option).",
                    path);
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

    void ZblPotential::loadDataset(std::istream& dataset, const std::set<Species2Sorted>* species_filter) {
        std::string line;
        while (std::getline(dataset, line)) {
            if (line.empty()) continue;
            std::istringstream iss(line);
            std::string elem1_str, elem2_str;
            std::array<Real, 6> coeffs = {};

            if (!(iss >> elem1_str >> elem2_str >> coeffs[0] >> coeffs[1] >> coeffs[2] >> coeffs[3] >> coeffs[4] >>
                  coeffs[5])) {
                continue;
            }

            Species s1(elem1_str);
            Species s2(elem2_str);
            Species2Sorted pair{s1, s2};

            if (species_filter && !species_filter->contains(pair)) {
                continue;
            }

            zbl_parameters[pair] = makeParameters(pair, coeffs);
        }
    }

    ZblPotential::ZblParameters ZblPotential::makeParameters(const Species2Sorted& pair,
                                                             const std::array<Real, 6>& coeffs) {
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

        return {coeffs, z1 * z2, (std::pow(z1, 0.23_r) + std::pow(z2, 0.23_r)) / 0.46848_r};
    }

    ZblPotential::ZblPotential(const std::set<Species2Sorted>& species, EmbeddedZBLCoeffDataset embedded_dataset,
                               Real cutoff, Real cutoff_transition_width) :
        cutoff(cutoff),
        cutoff_transition_width(cutoff_transition_width),
        cutoff_function(cutoff, cutoff_transition_width) {
        std::istringstream dataset_stream(getEmbeddedDataset(embedded_dataset));
        loadDataset(dataset_stream, &species);
    }

    ZblPotential::ZblPotential(const std::vector<Atoms>& training_data, EmbeddedZBLCoeffDataset embedded_dataset,
                               Real cutoff, Real cutoff_transition_width) :
        cutoff(cutoff),
        cutoff_transition_width(cutoff_transition_width),
        cutoff_function(cutoff, cutoff_transition_width) {
        std::set<Species2Sorted> pairs = deducePairs(training_data);
        std::istringstream dataset_stream(getEmbeddedDataset(embedded_dataset));
        loadDataset(dataset_stream, &pairs);
    }

    ZblPotential::ZblPotential(const std::string& dataset_filename, const std::vector<Atoms>& training_data,
                               Real cutoff, Real cutoff_transition_width) :
        cutoff(cutoff),
        cutoff_transition_width(cutoff_transition_width),
        cutoff_function(cutoff, cutoff_transition_width) {
        std::ifstream dataset_file(dataset_filename);
        if (!dataset_file.is_open()) {
            JGAP_LOG_AND_THROW("Could not open ZBL dataset file '{}'", dataset_filename);
        }
        std::set<Species2Sorted> pairs = deducePairs(training_data);
        loadDataset(dataset_file, &pairs);
    }

    ZblPotential::ZblPotential(std::istream& custom_dataset, const std::set<Species2Sorted>& species, Real cutoff,
                               Real cutoff_transition_width) :
        cutoff(cutoff),
        cutoff_transition_width(cutoff_transition_width),
        cutoff_function(cutoff, cutoff_transition_width) {
        loadDataset(custom_dataset, &species);
    }

    ZblPotential::ZblPotential(std::istream& custom_dataset, const std::vector<Atoms>& training_data, Real cutoff,
                               Real cutoff_transition_width) :
        cutoff(cutoff),
        cutoff_transition_width(cutoff_transition_width),
        cutoff_function(cutoff, cutoff_transition_width) {
        std::set<Species2Sorted> pairs = deducePairs(training_data);
        loadDataset(custom_dataset, &pairs);
    }

    ZblPotential::ZblPotential(const std::map<Species2Sorted, std::array<Real, 6>>& coefficients, Real cutoff,
                               Real cutoff_transition_width) :
        cutoff(cutoff),
        cutoff_transition_width(cutoff_transition_width),
        cutoff_function(cutoff, cutoff_transition_width) {
        for (const auto& [pair, coeffs]: coefficients) {
            zbl_parameters[pair] = makeParameters(pair, coeffs);
        }
    }

    std::map<Species2Sorted, std::array<Real, 6>> ZblPotential::getCoefficients() const {
        std::map<Species2Sorted, std::array<Real, 6>> coefficients;
        for (const auto& [pair, params]: zbl_parameters) {
            coefficients[pair] = params.coeffs;
        }
        return coefficients;
    }


    void ZblPotential::fillTables(TabulationData& tables) const {
        for (const auto& [species_pair, params]: zbl_parameters) {
            auto& table = tables.two_body_grids.getValueGrid(species_pair);
            for (const auto cell: table) {
                if (cell.pos[0] < LowerRLimitForTabulation) continue;
                cell.value += energyAndDerivative(params, cell.pos[0])[0];
            }
        }
    }

    AtomicQuantity ZblPotential::calculateEnergy(const Atoms& atoms) const {
        AtomicQuantity result(atoms.nAtoms());

        NeighbourLists nl(atoms, cutoff);

        for (const auto& [species_pair, params]: zbl_parameters) {
            Cluster2Expansion expansion(species_pair);
            auto expansion_result = expansion.find(nl, CalculationType::WithGradients);

            assert(expansion_result.derivatives.has_value());

            for (const auto& [cluster, cluster_derivs]:
                 std::views::zip(expansion_result.clusters, *expansion_result.derivatives)) {
                const Real& separation_magnitude = cluster.r01;
                const SeparationDerivatives& separation_deriv = cluster_derivs.dr01;

                auto [e, dE_dr] = energyAndDerivative(params, separation_magnitude);

                result.value += e;

                result.virials += separation_deriv.virials * dE_dr;

                Vector3 f01 = separation_deriv.direction * dE_dr;
                result.forces[cluster.atom_indexes[1]] -= f01;
                result.forces[cluster.atom_indexes[0]] += f01;
            }
        }

        return result;
    }

    std::array<Real, 2> ZblPotential::energyAndDerivative(const Species2Sorted& species_pair, Real r) const {
        auto it = zbl_parameters.find(species_pair);
        if (it == zbl_parameters.end()) {
            return {0.0, 0.0};
        }
        return energyAndDerivative(it->second, r);
    }

    std::array<Real, 2> ZblPotential::energyAndDerivative(const ZblParameters& params, Real r) const {
        Real x = r * params.a_inverse;
        Real dx_dr = params.a_inverse;

        Real term1 = exp(-params.coeffs[1] * x);
        Real term2 = exp(-params.coeffs[3] * x);
        Real term3 = exp(-params.coeffs[5] * x);

        Real phi = params.coeffs[0] * term1 + params.coeffs[2] * term2 + params.coeffs[4] * term3;

        Real dphi_dx = -params.coeffs[0] * params.coeffs[1] * term1 - params.coeffs[2] * params.coeffs[3] * term2 -
                       params.coeffs[4] * params.coeffs[5] * term3;

        Real dphi_dr = dx_dr * dphi_dx;

        Real prefactor = params.z1_z2 * CoulombConstant_eV_Ang / r;

        Real zbl_e = prefactor * phi;
        Real dzbl_dr = -(zbl_e / r) + prefactor * dphi_dr;

        auto [cut, dcut] = cutoff_function.evaluateAndDifferentiate(r);

        Real E = cut * zbl_e;
        Real dE_dr = cut * dzbl_dr + dcut * zbl_e;

        return {E, dE_dr};
    }
}
