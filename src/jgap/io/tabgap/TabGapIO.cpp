#include "jgap/io/tabgap/TabGapIO.hpp"

#include <algorithm>
#include <fstream>
#include <highfive/H5File.hpp>
#include <sstream>

#include "jgap/core/atomic/species/composition/Species2Atomic.hpp"
#include "jgap/core/potentials/tabgap/components/ThreeBodyTGComponent.hpp"
#include "jgap/core/potentials/tabgap/components/TwoBodyTGComponent.hpp"
#include "jgap/core/splines/CubicBSpline.hpp"
#include "jgap/core/transform/nbody/2b/eam/SplinePairTransformation.hpp"
#include "jgap/io/log/CurrentLogger.hpp"

namespace jgap {

    TabGapPotential TabGapIO::read(const Filenames& filenames) {
        size_t found_h5{};
        bool eam_fs_provided = false;
        for (auto& filename: filenames) {
            if (filename.ends_with(".h5") && found_h5++ > 1) {
                JGAP_LOG_AND_THROW(
                    "Multiple .h5 files detected when reading single tabulation data: {}", vectorToString(filenames)
                );
            }
            if (filename.ends_with(".eam.fs")) {
                eam_fs_provided = true;
            }
        }

        TabGapPotential potential{};
        for (const auto& filename: filenames) {
            JGAP_LOG_DEBUG("Reading file {}", filename);
            if (filename.ends_with(".h5")) {
                HighFive::File h5_file(filename, HighFive::File::ReadOnly);
                HighFive::Group root = h5_file.getGroup("/");

                // Fall back to the EAM data embedded in the .h5 if no .eam.fs files were supplied.
                readFromGroup(root, potential, !eam_fs_provided);
            } else if (filename.ends_with(".eam.fs")) {
                std::ifstream eam_fs_file(filename);
                if (!eam_fs_file.is_open()) {
                    JGAP_LOG_AND_THROW("Could not open {} as .eam.fs", filename);
                }
                parseEamFs(eam_fs_file, potential);
            } else {
                JGAP_LOG_AND_THROW("File format not supported to read as a tabGAP: {}", filename);
            }
        }

        potential.recomputeComponentCounts();
        return potential;
    }

    Filenames TabGapIO::write(const TabGapPotential& potential, const std::string& output_filename_prefix) {
        JGAP_LOG_INFO("Writing {} tabGAP files", output_filename_prefix);
        Filenames resulting_filenames{output_filename_prefix + ".tabgap.h5"};
        HighFive::File h5_file(resulting_filenames.back(), HighFive::File::Overwrite);
        HighFive::Group root = h5_file.getGroup("/");

        const std::vector<std::string> eam_fs_contents = writeToGroup(root, potential);
        h5_file.flush();

        if (eam_fs_contents.empty()) {
            return resulting_filenames;
        }

        JGAP_LOG_DEBUG("Writing {} .eam.fs file(s)", output_filename_prefix);
        for (size_t index{}; index < eam_fs_contents.size(); index++) {
            const std::string filename =
                output_filename_prefix + (index != 0 ? "#" + std::to_string(index) : "") + ".eam.fs";

            std::ofstream eam_fs_file(filename);
            if (!eam_fs_file.is_open()) {
                JGAP_LOG_ERROR("Could not open {}, saving to H5 only", filename);
                continue;
            }
            eam_fs_file << eam_fs_contents[index];
            eam_fs_file.flush();

            resulting_filenames.push_back(filename);
        }

        return resulting_filenames;
    }

    std::vector<std::string> TabGapIO::writeToGroup(HighFive::Group& root, const TabGapPotential& potential) {
        const std::string comment1 = "UNITS: metal";
        root.createDataSet<std::string>("comment1", comment1);
        const std::string comment2 = "pair_style tabgap";
        root.createDataSet<std::string>("comment2", comment2);

        auto e0Group = root.createGroup("e0");
        e0Group.createAttribute("Nelements", potential.isolated_atom_energies.size());
        for (const auto& [element, value]: potential.isolated_atom_energies) {
            e0Group.createAttribute(element.symbol(), value);
        }

        if (potential.n_eam_components == 0) {
            root.createDataSet("npots", std::array{potential.n_2b_components, potential.n_3b_components});
        } else {
            root.createDataSet("npots", std::array{size_t{0}, potential.n_3b_components});
        }

        std::set<Species> species_2b_and_eam;
        std::map<Species2Sorted, const TwoBodyTGComponent*> pair_pots;
        std::multimap<Species, const EamTGComponent*> eam_components;
        for (auto& component: potential.components) {
            if (const auto* casted2b = dynamic_cast<const TwoBodyTGComponent*>(component.get()); casted2b != nullptr) {
                if (potential.n_eam_components == 0) {
                    write2b(root, *casted2b);
                } else {
                    pair_pots.insert({casted2b->getSpeciesPair(), casted2b});
                    species_2b_and_eam.insert(casted2b->getSpeciesPair().nodes[0]);
                    species_2b_and_eam.insert(casted2b->getSpeciesPair().nodes[1]);
                }

            } else if (
                const auto* casted3b = dynamic_cast<const ThreeBodyTGComponent*>(component.get()); casted3b != nullptr
            ) {
                write3b(root, *casted3b);

            } else if (
                const auto* casted_eam = dynamic_cast<const EamTGComponent*>(component.get()); casted_eam != nullptr
            ) {
                eam_components.insert({casted_eam->getSplineNBodyAggregator().getCentralSpecies(), casted_eam});

                for (Species s: casted_eam->getSplineNBodyAggregator().getAllSpecies()) {
                    species_2b_and_eam.insert(s);
                }

            } else {
                JGAP_LOG_AND_THROW("tabGAP component of type not supported for serialization");
            }
        }

        std::vector<std::string> eam_fs_contents;
        if (potential.n_eam_components) {
            std::vector species_2b_and_eam_vec(species_2b_and_eam.begin(), species_2b_and_eam.end());

            while (!pair_pots.empty() || !eam_components.empty()) {
                eam_fs_contents.push_back(
                    useSomeComponentsAndGenerateEamFs(species_2b_and_eam_vec, pair_pots, eam_components)
                );
            }

            // Embed the .eam.fs contents as root attributes so the potential is self-contained in the .h5
            // without adding top-level groups that third-party tabgap Python tools iterate over.
            root.createAttribute("Neam_files", eam_fs_contents.size());
            for (size_t index{}; index < eam_fs_contents.size(); index++) {
                root.createAttribute(std::format("eam_file_{}", index), eam_fs_contents[index]);
            }
        }

        return eam_fs_contents;
    }

    /// Writes one 2-body group into `root`.
    ///
    /// On-disk convention (a quirk of the external tabGAP format): the group describes the *original*
    /// tabulation grid but stores the cubic B-spline *coefficients* as its values:
    /// <ul>
    ///   <li>`N`: original grid point count (= coefficient count - 2; CubicBSpline::fit pads one
    ///       ghost coefficient on each side).</li>
    ///   <li>`grid_limits`: [lower, upper]. `lower` is the original-grid origin
    ///       (coeff origin + spacing); `upper` is the coefficient grid's cutoff, which sits ONE
    ///       spacing past the last original point (origin + (coeff_dims - 1) * spacing).</li>
    ///   <li>`energies`: the raw B-spline coefficients (N + 2 values).</li>
    /// </ul>
    /// Because of the one-extra-spacing `upper`, the reader recovers spacing as (upper - lower) / N
    /// (see \ref TabGapIO::readFromGroup), not / (N - 1).
    ///
    /// @param root the HDF5 group to add the pair group to.
    /// @param component the 2-body component to write (its spline must be a CubicBSpline).
    void TabGapIO::write2b(HighFive::Group& root, const TwoBodyTGComponent& component) {
        auto species = component.getSpeciesPair().nodes;

        auto pair_group = root.createGroup(std::format("{}-{}", species[0].symbol(), species[1].symbol()));

        pair_group.createAttribute("element_i", species[0].symbol());
        pair_group.createAttribute("element_j", species[1].symbol());

        auto spline_cast = dynamic_cast<const CubicBSpline*>(component.getSpline().get());
        if (!spline_cast) {
            JGAP_LOG_AND_THROW("2b tabGAP component must use CubicBSpline to be writeable into H5");
        }

        const auto& coeff_grid = spline_cast->getCoefficients();
        const auto original_grid_size = static_cast<Real>(coeff_grid.sizes[0] - 2);

        // Original grid bounds: lower = origin + spacing, upper = origin + (sizes - 2) * spacing
        pair_group.createDataSet(/*Original grid, not spline coeffs*/
                                 "grid_limits",
                                 std::vector{
                                     coeff_grid.origin[0] + coeff_grid.spacing[0],
                                     coeff_grid.origin[0] + original_grid_size * coeff_grid.spacing[0]
                                 });

        pair_group.createAttribute("N", coeff_grid.sizes[0] - 2 /*Original grid, not spline coeffs*/);
        pair_group.createDataSet("energies", coeff_grid.data_flat);
    }

    /// Writes one 3-body group, following the same on-disk convention as \ref TabGapIO::write2b (see its doc):
    /// per axis, `N` is the original grid point count (= coeff dims - 2) and `grid_limits`
    /// stores the original-grid extent [origin, cutoff] derived from the coefficient grid; `energies`
    /// holds the raw B-spline coefficients (the original count + 2 per axis). Each triplet group is named
    /// "<i>-<j>-<k>" so readers can tell 2b (one '-') from 3b (two '-') groups.
    ///
    /// @param root the HDF5 group to add the triplet group to.
    /// @param component the 3-body component to write.
    void TabGapIO::write3b(HighFive::Group& root, const ThreeBodyTGComponent& component) {
        auto [root_species, node_species] = component.getSpeciesTriplet();

        auto triplet_group = root.createGroup(
            std::format("{}-{}-{}", root_species.symbol(), node_species[0].symbol(), node_species[1].symbol())
        );

        triplet_group.createAttribute("element_i", root_species.symbol());
        triplet_group.createAttribute("element_j", node_species[0].symbol());
        triplet_group.createAttribute("element_k", node_species[1].symbol());

        const auto& coeff_grid = component.getSpline().getCoefficients();
        const auto original_grid_sizes = std::array{
            static_cast<Real>(coeff_grid.sizes[0] - 2),
            static_cast<Real>(coeff_grid.sizes[1] - 2),
            static_cast<Real>(coeff_grid.sizes[2] - 2),
        };

        // Original-grid lower limit = coeff origin + spacing; upper limit = coeff origin + (sizes - 2) * spacing
        triplet_group.createDataSet(
            "grid_limits", std::array{
                               coeff_grid.origin[0] + coeff_grid.spacing[0],
                               coeff_grid.origin[1] + coeff_grid.spacing[1],
                               coeff_grid.origin[2] + coeff_grid.spacing[2],
                               coeff_grid.origin[0] + original_grid_sizes[0] * coeff_grid.spacing[0],
                               coeff_grid.origin[1] + original_grid_sizes[1] * coeff_grid.spacing[1],
                               coeff_grid.origin[2] + original_grid_sizes[2] * coeff_grid.spacing[2],
                           }
        );

        triplet_group.createDataSet(
            "N", std::array{
                     // original grid point counts N, not coeff counts (coeff count = N + 2)
                     coeff_grid.sizes[0] - 2, coeff_grid.sizes[1] - 2, coeff_grid.sizes[2] - 2
                 }
        );

        triplet_group.createDataSet("energies", coeff_grid.data_flat);
    }

    /// Builds the contents of one LAMMPS .eam.fs file from (some of) the supplied components, and
    /// CONSUMES the ones it used: it takes at most one EAM component per element plus the pair potentials,
    /// erasing them from `eam_components`/`pair_pots`. Callers loop until both are empty,
    /// producing one .eam.fs per round (a single .eam.fs cannot hold two embedding functions for the same
    /// element).
    ///
    /// .eam.fs structure: 3 comment lines; an elements line; a (Nrho, drho, Nr, dr, cutoff) line; then per
    /// element an embedding function F(rho) and one density function rho(r) per element; finally the
    /// pair-potential tables phi(r) for i &gt;= j (stored as r * phi(r), as the format requires).
    ///
    /// @param all_species the element ordering shared by every emitted file.
    /// @param pair_pots pair potentials still to write; consumed (cleared) by this call.
    /// @param eam_components EAM components still to write; the ones used are erased.
    /// @return the text of one .eam.fs file.
    std::string TabGapIO::useSomeComponentsAndGenerateEamFs(
        const std::vector<Species>& all_species, std::map<Species2Sorted, const TwoBodyTGComponent*>& pair_pots,
        std::multimap<Species, const EamTGComponent*>& eam_components
    ) {
        ////////// Pre-process ///////////
        if (eam_components.empty()) {
            JGAP_LOG_AND_THROW("Unexpected behaviour while writing .eam.fs");
        }

        size_t n_rho{}, n_2b{};
        Real drho{}, dr{};
        std::map<Species, Grid<1>> energy_per_density_grids;
        std::map<Species2Atomic, Grid<1>> density_grids;

        for (Species element: all_species) {
            if (eam_components.contains(element)) {
                auto it = eam_components.find(element);

                energy_per_density_grids.insert({element, it->second->getEnergySpline().getTable()});
                drho = it->second->getEnergySpline().getTable().spacing[0];
                n_rho = it->second->getEnergySpline().getTable().sizes[0];

                for (const auto& [species_pair, eam_pf_trans]:
                     it->second->getSplineNBodyAggregator().getTransformations()) {
                    if (auto as_spline_trans = dynamic_cast<const SplinePairTransformation*>(eam_pf_trans.get())) {
                        density_grids.insert({species_pair, as_spline_trans->getSpline().getTable()});

                        n_2b = as_spline_trans->getSpline().getTable().sizes[0];
                        dr = as_spline_trans->getSpline().getTable().spacing[0];
                    } else {
                        JGAP_LOG_AND_THROW("Non spline aggregator detected");
                    }
                }

                eam_components.erase(it);
            }
        }

        // not a ref:
        std::map<Species2Sorted, Grid<1>> pair_pot_grids{};

        for (const auto& [species_pair, pair_pot]: pair_pots) {
            if (auto as_hermite = dynamic_cast<const HermiteCubicSpline*>(pair_pot->getSpline().get());
                as_hermite != nullptr) {
                pair_pot_grids.insert({species_pair, as_hermite->getTable()});

            } else {
                JGAP_LOG_AND_THROW("Pairpot spline must be HermiteCubic to be saved in eam.fs");
            }
        }

        pair_pots.clear();

        ////////// Actually write ///////////

        std::ostringstream eam_fs_content;
        eam_fs_content << std::fixed << std::setprecision(17);

        // Lines 1–3: Comments/metadata.
        eam_fs_content << "# UNITS: metal" << std::endl;
        eam_fs_content << "# EAM part of a potential tabulated with jGAP" << std::endl;
        eam_fs_content << "# pair_style eam/fs" << std::endl;

        // Line 4: Number of elements (N) followed by each element’s symbol
        eam_fs_content << all_species.size() << " ";
        for (const auto& element: all_species) {
            eam_fs_content << element.symbol() << " ";
        }
        eam_fs_content << std::endl;

        // Line 5: Nrho, drho, Nr, dr, cutoff - so mysterious ..
        eam_fs_content << n_rho << " ";
        eam_fs_content << drho << " ";
        eam_fs_content << n_2b << " ";
        eam_fs_content << dr << " ";
        eam_fs_content << dr * static_cast<Real>(n_2b) << std::endl;

        // Per-element Sections:
        /*
         * Line: atomic number, mass, lattice constant, lattice type (e.g., fcc, bcc)
         * Embedding function F_\beta(\rho): Nrho tabulated values
         * Density functions \rho_{\alpha\beta}(r): For each element α (total N curves, each with Nr points)
         */
        for (const Species& contributor: all_species) {
            eam_fs_content << contributor.atomicNumber().value_or(-1) << " ";
            eam_fs_content << contributor.mass().value_or(-1) << " 1.0 ZZZ" << std::endl;

            if (energy_per_density_grids.contains(contributor)) {
                for (const auto& density_grid_cell: energy_per_density_grids.at(contributor)) {
                    eam_fs_content << density_grid_cell.value << std::endl;
                }
            } else {
                for (size_t i = 0; i < n_rho; i++) {
                    eam_fs_content << 0.0 << std::endl;
                }
            }

            // contributor = species1
            for (const Species& receiver: all_species) {
                if (density_grids.contains(Species2Atomic{receiver, contributor})) {
                    for (const auto& density_grid_cell: density_grids.at(Species2Atomic{receiver, contributor})) {
                        eam_fs_content << density_grid_cell.value << std::endl;
                    }
                } else {
                    for (size_t i = 0; i < n_2b; i++) {
                        eam_fs_content << 0.0 << std::endl;
                    }
                }
            }
        }

        /*
         * Pair Potential Tables (for all i ≥ j pairs):
         * Tabulated \phi_{\alpha\beta}(r) values for each unique pair (symmetry is exploited),
         * listing Nr points per interaction
         */

        for (size_t i = 0; i < all_species.size(); i++) {
            for (size_t j = 0; j < all_species.size(); j++) {
                if (i < j) {
                    continue;
                }

                Species2Sorted species_pair{all_species[i], all_species[j]};

                if (pair_pot_grids.contains(species_pair)) {
                    for (auto cell: pair_pot_grids.at(species_pair)) {
                        eam_fs_content << cell.value * cell.pos[0] << std::endl;
                    }
                } else {
                    for (size_t k = 0; k < n_2b; k++) {
                        eam_fs_content << 0.0 << std::endl;
                    }
                }
            }
        }

        return eam_fs_content.str();
    }

    TabGapPotential TabGapIO::fromGroup(const HighFive::Group& root, bool read_embedded_eam_fs) {
        TabGapPotential pot{};
        readFromGroup(root, pot, read_embedded_eam_fs);
        pot.recomputeComponentCounts();
        return pot;
    }

    void TabGapIO::readFromGroup(const HighFive::Group& root, TabGapPotential& pot, bool read_embedded_eam_fs) {
        // Read isolated atom energies
        if (root.exist("e0")) {
            auto e0_group = root.getGroup("e0");
            // Read all attributes except Nelements
            for (const auto& attr_name: e0_group.listAttributeNames()) {
                if (attr_name == "Nelements") continue;
                Real val;
                e0_group.getAttribute(attr_name).read(val);
                pot.isolated_atom_energies[attr_name] += val;
            }
        }

        // Iterate groups by name to find 2b/3b components.
        for (const auto& name: root.listObjectNames()) {
            if (name == "e0" || name == "npots" || name == "comment1" || name == "comment2" || name == "eam_files") {
                continue;
            }
            // Count number of '-' to distinguish 2b vs 3b
            size_t dash_count = static_cast<size_t>(std::ranges::count(name, '-'));
            auto group = root.getGroup(name);

            if (dash_count == 1) {
                // Pair group: element_i, element_j
                std::string species_i, species_j;
                group.getAttribute("element_i").read(species_i);
                group.getAttribute("element_j").read(species_j);
                Species2Sorted pair{species_i, species_j};

                std::array<Real, 2> limits{}; // origin, cutoff
                group.getDataSet("grid_limits").read(limits);

                size_t n_original_grid;
                group.getAttribute("N").read(n_original_grid);

                // Inverse of write2b: limits store original grid bounds [lower, upper] = [r_min, r_max].
                // Original grid has N points spanning [lower, upper], so spacing = (upper - lower) / (N - 1).
                // Spline grid has N + 2 coefficient points with origin shifted back by 1 spacing (lower - spacing).
                Real lower = limits.at(0);
                Real upper = limits.at(1);
                Real spacing = (upper - lower) / static_cast<Real>(n_original_grid - 1);
                Grid<1> spline_grid{{n_original_grid + 2}, {spacing}, {lower - spacing}};

                group.getDataSet("energies").read(spline_grid.data_flat);
                if (spline_grid.data_flat.size() != n_original_grid + 2) {
                    JGAP_LOG_AND_THROW("Invalid H5 pair data size for {}-{}", species_i, species_j);
                }

                pot.components.emplace_back(TwoBodyTGComponent{pair, CubicBSpline(spline_grid)});

            } else if (dash_count == 2) {
                // Triplet group: element_i, element_j, element_k
                std::string species_i, species_j, species_k;
                group.getAttribute("element_i").read(species_i);
                group.getAttribute("element_j").read(species_j);
                group.getAttribute("element_k").read(species_k);

                Species3AtomicSorted triplet{species_i, species_j, species_k};

                std::array<size_t, 3> n_original{}; // original grid point counts N per axis (see write3b)
                group.getDataSet("N").read(n_original);

                std::array<Real, 6> grid_limits{}; // lower xyz, then upper xyz
                group.getDataSet("grid_limits").read(grid_limits);

                // Inverse of write3b: limits store original grid bounds [lower, upper] per axis.
                // Original grid has N points spanning [lower, upper], so spacing = (upper - lower) / (N - 1).
                // Spline grid has N + 2 coefficient points with origin shifted back by 1 spacing (lower - spacing).
                std::array lower{grid_limits[0], grid_limits[1], grid_limits[2]};
                std::array upper{grid_limits[3], grid_limits[4], grid_limits[5]};
                std::array spacing{
                    (upper[0] - lower[0]) / static_cast<Real>(n_original[0] - 1),
                    (upper[1] - lower[1]) / static_cast<Real>(n_original[1] - 1),
                    (upper[2] - lower[2]) / static_cast<Real>(n_original[2] - 1)
                };
                std::array<size_t, 3> spline_dims{n_original[0] + 2, n_original[1] + 2, n_original[2] + 2};
                std::array spline_grid_origin{lower[0] - spacing[0], lower[1] - spacing[1], lower[2] - spacing[2]};

                std::vector<Real> spline_coeffs{};
                group.getDataSet("energies").read(spline_coeffs);

                assert(
                    spline_coeffs.size() == spline_dims[0] * spline_dims[1] * spline_dims[2] &&
                    "Number of coefficients mis-matches grid size"
                );

                Grid spline_grid{spline_dims, spacing, spline_grid_origin, spline_coeffs};

                pot.components.emplace_back(ThreeBodyTGComponent(triplet, CubicBSpline3D(spline_grid)));
            }
        }

        // EAM part: take it from embedded root attributes when asked to.
        if (read_embedded_eam_fs) {
            if (root.hasAttribute("Neam_files")) {
                size_t n_eam_files{};
                root.getAttribute("Neam_files").read(n_eam_files);
                for (size_t index = 0; index < n_eam_files; index++) {
                    std::string eam_fs_content;
                    root.getAttribute(std::format("eam_file_{}", index)).read(eam_fs_content);
                    std::istringstream eam_fs_stream(eam_fs_content);
                    parseEamFs(eam_fs_stream, pot);
                }
            } else if (root.exist("eam_files")) {
                auto eam_files_group = root.getGroup("eam_files");
                for (const auto& eam_fs_name: eam_files_group.listObjectNames()) {
                    std::string eam_fs_content;
                    eam_files_group.getDataSet(eam_fs_name).read(eam_fs_content);
                    std::istringstream eam_fs_stream(eam_fs_content);
                    parseEamFs(eam_fs_stream, pot);
                }
            }
        }
    }

    /// Parses one LAMMPS .eam.fs stream (a file or an embedded "eam_files" dataset) and appends the
    /// resulting EAM components (one per element) and pair-potential 2-body components to `pot`.
    /// Inverse of \ref ::useSomeComponentsAndGenerateEamFs; pair tables stored as r * phi(r) are divided
    /// back by r.
    ///
    /// @param file the .eam.fs text stream to parse.
    /// @param pot the potential to append the parsed components to.
    void TabGapIO::parseEamFs(std::istream& file, TabGapPotential& pot) {
        std::string line;
        for (size_t i = 0; i < 3; i++) {
            if (!getLine(file, line)) JGAP_LOG_AND_THROW("Invalid EAM/FS: missing comment #{}", i);
        }

        // Elements line: N and symbols
        if (!getLine(file, line)) JGAP_LOG_AND_THROW("Invalid EAM/FS: missing elements line");

        std::istringstream iss(line);
        size_t N;
        iss >> N;
        std::vector<Species> elements;
        for (size_t i = 0; i < N; i++) {
            std::string element_str;
            iss >> element_str;

            if (std::ranges::find(elements, Species(element_str)) != elements.end()) {
                JGAP_LOG_AND_THROW("Invalid EAM/FS: duplicate element {}", element_str);
            }

            elements.emplace_back(element_str);
        }
        if (elements.empty()) JGAP_LOG_AND_THROW("Invalid EAM/FS: zero elements");

        // Grid sizes line: Nrho, drho, Nr, dr, cutoff
        if (!getLine(file, line)) JGAP_LOG_AND_THROW("Invalid EAM/FS: missing grid spec line");

        size_t n_rho, n_r;
        Real drho, dr, cutoff;
        iss = std::istringstream(line);
        iss >> n_rho >> drho >> n_r >> dr >> cutoff;

        std::map<Species, Grid<1>> embedding_energies;
        std::map<Species, AtomicTwoBodyGrids<1>> density_grids;

        for (const auto& central_atom_species: elements) {
            density_grids.insert({central_atom_species, AtomicTwoBodyGrids<1>{{0.0}, {dr}, {n_r}}});
        }

        // Per-element sections
        for (const auto& central_atom_species: elements) {
            // Header line with atomic number, mass, etc. Ignore contents
            if (!getLine(file, line)) {
                JGAP_LOG_AND_THROW("Invalid EAM/FS: missing per-element header");
            }

            Grid<1> energies{{n_rho}, {drho}, {0.0}};
            // Embedding function: Nrho lines
            for (size_t i = 0; i < n_rho; i++) {
                if (!getLine(file, line)) JGAP_LOG_AND_THROW("Invalid EAM/FS: incomplete embedding values");
                energies.data_flat[i] = std::stod(line);
            }
            embedding_energies.insert({central_atom_species, std::move(energies)});

            // eam/fs convention (LAMMPS): in element β's section the N listed densities are rho_{αβ}(r),
            // the density contributed BY a β neighbour at an α atom (β = central_atom_species here, the
            // section element; α = receiver_species, the inner index). Such a density belongs to the EAM
            // component for CENTRAL atom α with NEIGHBOUR β, i.e. Species2Sorted{α, β}.
            for (const auto& receiver_species: elements) {
                Grid<1> rho_per_r_grid{{n_r}, {dr}, {0.0}};

                for (size_t i = 0; i < n_r; i++) {
                    if (!getLine(file, line)) {
                        JGAP_LOG_AND_THROW("Invalid EAM/FS: incomplete density function");
                    }
                    rho_per_r_grid.data_flat[i] = std::stod(line);
                }

                density_grids.at(receiver_species).getValueGrid({receiver_species, central_atom_species}) =
                    std::move(rho_per_r_grid);
            }
        }

        for (auto central_species: elements) {
            pot.components.emplace_back(EamTGComponent(
                ManyBodyGrids2<1, 1>{
                    .central_atom_species = central_species,
                    .aggregator_grids = std::move(density_grids.at(central_species)),
                    .value_grid = std::move(embedding_energies[central_species])
                }
            ));
        }

        // Pair potential tables for i >= j
        for (size_t i = 0; i < N; i++) {
            for (size_t j = 0; j < N; j++) {
                if (i < j) continue;

                Species2Sorted species_pair{elements[i], elements[j]};
                Grid<1> energy_grid({n_r}, {dr}, {0.0});

                Real r{};
                for (size_t k = 0; k < n_r; k++, r += dr) {
                    if (!getLine(file, line)) {
                        JGAP_LOG_AND_THROW("Invalid EAM/FS: incomplete pair potential table");
                    }
                    Real phi = std::stod(line);
                    energy_grid.data_flat[k] = (r > 0.0 ? phi / r : 0.0);
                }

                pot.components.emplace_back(TwoBodyTGComponent{species_pair, HermiteCubicSpline{energy_grid}});
            }
        }
    }
}
