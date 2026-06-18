#include "ZblPotential.hpp"
#include "core/cutoff/PerriotPolynomialCutoff.hpp"

#include <fstream>
#include <filesystem>
#include <sstream>
#include <cmath>

#include "core/atomic/neighbours/NeighbourList.hpp"
#include "ZblResources.h"

namespace jgap {

    namespace {
#if defined(__has_embed)
        std::string getEmbeddedDataset(EmbeddedZBLCoeffDataset type) {
            switch (type) {
                case EmbeddedZBLCoeffDataset::NLH:
                    return std::string(reinterpret_cast<const char*>(nlh_dat), nlh_dat_len);
                case EmbeddedZBLCoeffDataset::DMOL:
                    return std::string(reinterpret_cast<const char*>(dmol_dat), dmol_dat_len);
                case EmbeddedZBLCoeffDataset::MP2:
                    return std::string(reinterpret_cast<const char*>(mp2_dat), mp2_dat_len);
            }
            return "";
        }
#endif

        std::set<SpeciesSet<2, Symmetric>> deducePairs(const std::vector<Atoms>& training_data) {
            std::set<SpeciesSet<2, Symmetric>> species_pairs;
            for (const auto& atoms : training_data) {
                const auto& species_vec = atoms.lookupSpecies();
                std::set unique_species(species_vec.begin(), species_vec.end());
                for (auto it1 = unique_species.begin(); it1 != unique_species.end(); ++it1) {
                    for (auto it2 = it1; it2 != unique_species.end(); ++it2) {
                        species_pairs.insert(SpeciesSet<2, Symmetric>{*it1, *it2});
                    }
                }
            }
            return species_pairs;
        }
    }

    void ZblPotential::loadDataset(std::istream& dataset, const std::set<SpeciesSet<2, Symmetric>>* species_filter) {
        std::string line;
        while (std::getline(dataset, line)) {
            if (line.empty()) continue;
            std::istringstream iss(line);
            std::string elem1_str, elem2_str;
            std::array<Real, 6> coeffs = {};

            if (!(iss >> elem1_str >> elem2_str >> coeffs[0] >> coeffs[1] >> coeffs[2] >> coeffs[3] >> coeffs[4] >> coeffs[5])) {
                continue;
            }

            Species s1(elem1_str);
            Species s2(elem2_str);
            SpeciesSet<2, Symmetric> pair{s1, s2};

            if (species_filter && !species_filter->contains(pair)) {
                continue;
            }

            zbl_parameters[pair] = makeParameters(pair, coeffs);
        }
    }

    ZblPotential::ZblParameters ZblPotential::makeParameters(const SpeciesSet<2, Symmetric>& pair,
                                                             const std::array<Real, 6>& coeffs) {
        const Species& s1 = pair.getNodes()[0];
        const Species& s2 = pair.getNodes()[1];

        if (!s1.atomicNumber().has_value()) {
            JGAP_LOG_AND_THROW("Atomic number undefined for {}", s1.symbol());
        }
        if (!s2.atomicNumber().has_value()) {
            JGAP_LOG_AND_THROW("Atomic number undefined for {}", s2.symbol());
        }

        Real z1 = s1.atomicNumber().value();
        Real z2 = s2.atomicNumber().value();

        return {
            coeffs,
            z1 * z2,
            (std::pow(z1, 0.23) + std::pow(z2, 0.23)) / 0.46848
        };
    }

#if defined(__has_embed)
    ZblPotential::ZblPotential(const std::set<SpeciesSet<2, Symmetric>>& species,
                               EmbeddedZBLCoeffDataset embedded_dataset,
                               Real cutoff,
                               Real cutoff_transition_width)
        : cutoff(cutoff), cutoff_transition_width(cutoff_transition_width),
          cutoff_function(cutoff, cutoff_transition_width) {
        std::istringstream dataset_stream(getEmbeddedDataset(embedded_dataset));
        loadDataset(dataset_stream, &species);
    }

    ZblPotential::ZblPotential(const std::vector<Atoms>& training_data,
                               EmbeddedZBLCoeffDataset embedded_dataset,
                               Real cutoff,
                               Real cutoff_transition_width)
        : cutoff(cutoff), cutoff_transition_width(cutoff_transition_width),
          cutoff_function(cutoff, cutoff_transition_width) {
        std::set<SpeciesSet<2, Symmetric>> pairs = deducePairs(training_data);
        std::istringstream dataset_stream(getEmbeddedDataset(embedded_dataset));
        loadDataset(dataset_stream, &pairs);
    }
#endif

    ZblPotential::ZblPotential(std::istream& custom_dataset,
                               const std::set<SpeciesSet<2, Symmetric>>& species,
                               Real cutoff,
                               Real cutoff_transition_width)
        : cutoff(cutoff), cutoff_transition_width(cutoff_transition_width),
          cutoff_function(cutoff, cutoff_transition_width) {
        loadDataset(custom_dataset, &species);
    }

    ZblPotential::ZblPotential(std::istream& custom_dataset,
                               const std::vector<Atoms>& training_data,
                               Real cutoff,
                               Real cutoff_transition_width)
        : cutoff(cutoff), cutoff_transition_width(cutoff_transition_width),
          cutoff_function(cutoff, cutoff_transition_width) {
        std::set<SpeciesSet<2, Symmetric>> pairs = deducePairs(training_data);
        loadDataset(custom_dataset, &pairs);
    }

    ZblPotential::ZblPotential(const std::map<SpeciesSet<2, Symmetric>, std::array<Real, 6>>& coefficients,
                               Real cutoff,
                               Real cutoff_transition_width)
        : cutoff(cutoff), cutoff_transition_width(cutoff_transition_width),
          cutoff_function(cutoff, cutoff_transition_width) {
        for (const auto& [pair, coeffs] : coefficients) {
            zbl_parameters[pair] = makeParameters(pair, coeffs);
        }
    }

    std::map<SpeciesSet<2, Symmetric>, std::array<Real, 6>> ZblPotential::getCoefficients() const {
        std::map<SpeciesSet<2, Symmetric>, std::array<Real, 6>> coefficients;
        for (const auto& [pair, params] : zbl_parameters) {
            coefficients[pair] = params.coeffs;
        }
        return coefficients;
    }


    void ZblPotential::fillTables(TabulationData &tables) const {
        for (const auto& [species_pair, params]: zbl_parameters) {
            auto& table = tables.two_body_grids.getValueGrid(species_pair);
            for (const auto cell: table) {
                if (cell.pos[0] < LowerRLimitForTabulation) continue;
                cell.value += energyAndDerivative(params, cell.pos[0])[0];
            }
        }
    }

    AtomicQuantity ZblPotential::calculateEnergy(const Atoms &atoms) const {
        AtomicQuantity result(atoms.nAtoms());

        NeighbourList nl(atoms, cutoff);

        for (const auto &[species_pair, params]: zbl_parameters) {
            auto pair_clusters = nl.findAllClusters<WithGradients>(species_pair);

            for (auto pair: pair_clusters) {
                const Real& separation_magnitude = pair.between(0, 1);
                const SeparationDerivatives& separation_deriv = pair.derivativesBetween(0, 1);

                auto [e, dE_dr] = energyAndDerivative(params, separation_magnitude);

                result.value += e;

                result.virials += separation_deriv.virials * dE_dr;

                Vector3 f01 = separation_deriv.direction * dE_dr;
                result.forces[pair.atom_indexes[1]] -= f01;
                result.forces[pair.atom_indexes[0]] += f01;
            }
        }

        return result;
    }

    std::array<Real, 2> ZblPotential::energyAndDerivative(const SpeciesSet<2, Symmetric>& species_pair, Real r) const {
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

        Real phi = params.coeffs[0] * term1
                 + params.coeffs[2] * term2
                 + params.coeffs[4] * term3;

        Real dphi_dx = - params.coeffs[0] * params.coeffs[1] * term1
                       - params.coeffs[2] * params.coeffs[3] * term2
                       - params.coeffs[4] * params.coeffs[5] * term3;

        Real dphi_dr = dx_dr * dphi_dx;

        Real prefactor = params.z1_z2 * CoulombConstant_eV_Ang / r;

        Real zbl_e = prefactor * phi;
        Real dzbl_dr = - (zbl_e / r) + prefactor * dphi_dr;

        auto [cut, dcut] = cutoff_function.evaluateAndDifferentiate(r);

        Real E = cut * zbl_e;
        Real dE_dr = cut * dzbl_dr + dcut * zbl_e;

        return {E, dE_dr};
    }
}