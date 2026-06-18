#ifndef JGAP_QUIPXMLCONVERTER_HPP
#define JGAP_QUIPXMLCONVERTER_HPP

#include <pugixml.hpp>

#include "core/potentials/Potential.hpp"
#include "core/potentials/gap/GapPotential.hpp"
#include "core/potentials/isolated/IsolatedAtomPotential.hpp"

namespace jgap {
    /**
     * Converts a QUIP potential, parsed from its `.xml`, into a jgap {@link Potential}.
     *
     * Handles the GAP part (a {@code GAP_params} block: isolated-atom e0 references plus the
     * {@code gpSparse} sparse descriptors for 2-body, 3-body and EAM terms) and the "IP Glue" pair
     * potential ({@code pairpot}). The result mirrors how the terms compose: a {@link GapPotential} whose
     * external potential carries the isolated-atom energies and (if present) the glue pair potential,
     * possibly wrapped in a {@link CompositePotential}.
     *
     * This depends on pugixml and is therefore only compiled when XML support is enabled (the {@code xml}
     * vcpkg feature; see the README and CMakeLists). It reads the external {@code sparseX} files that the
     * .xml references by name, resolved relative to the current working directory.
     *
     * Note: tightly coupled to the QUIP .xml layout; format changes will surface as parse errors or
     * exceptions rather than silent misreads.
     */
    class QuipXmlConverter {
    public:
        /**
         * Converts a QUIP potential node into a {@link Potential}. Accepts a {@code GAP_params} node, a
         * {@code pairpot} node, or a wrapper node containing either/both as children; when both are
         * present they are combined (the glue pair potential becomes the GAP's external potential, via a
         * {@link CompositePotential} if an isolated-atom potential is also present).
         *
         * @param quip_potential_encoded the parsed XML node (e.g. the document root element).
         * @return the converted potential.
         * @throws std::runtime_error if no recognised potential is found in the node.
         */
        static ValuePtr<Potential> transform(const pugi::xml_node& quip_potential_encoded);

        /**
         * Converts an "IP Glue" {@code pairpot} node (a {@code Glue_params} block of per-pair (r, E)
         * point tables) into a {@link SplinePairPotential}.
         *
         * @param quip_pairpot the {@code pairpot} node.
         * @return the pair potential.
         */
        static ValuePtr<Potential> transformPairpot(const pugi::xml_node& quip_pairpot);

        /**
         * Converts a {@code GAP_params} node into a {@link GapPotential}: the {@code GAP_data} e0 values
         * become an {@link IsolatedAtomPotential} (set as the GAP's external potential) and the
         * {@code gpSparse} descriptors become the GAP components. If there is no {@code gpSparse}, only
         * the isolated-atom potential is returned.
         *
         * @param quip_gap_params the {@code GAP_params} node.
         * @return the GAP potential (or just the isolated-atom potential).
         */
        static ValuePtr<Potential> transformGapParams(const pugi::xml_node& quip_gap_params);

    private:
        /** Descriptor header fields parsed from a {@code gpCoordinates} descriptor string. */
        struct QuipDescriptorData {
            std::string type;
            double delta;
            double theta;
            double cutoff;

            std::optional<double> r_min;
            std::optional<double> cutoff_transition_width;
            std::optional<std::string> pair_function;
            std::optional<double> order;
            std::optional<std::string> mode;
        };

        /** Reads the {@code GAP_data} e0 entries into an {@link IsolatedAtomPotential}. */
        static IsolatedAtomPotential transformIsolatedAtomParams(
                                                        const pugi::xml_node& quip_isolated_atom_params);

        /** Builds the GAP components from the {@code gpSparse} block (one per {@code gpCoordinates}). */
        static GapPotential transformSparseData(const pugi::xml_node& quip_sparse_data);

        /**
         * Builds a 2-body component (TwoBodyTransformation + SquaredExpKernel) from a {@code distance_2b}
         * descriptor, reading sparse points from the referenced {@code sparseX} file and the per-point
         * {@code alpha} coefficients.
         */
        static ValuePtr<GapComponent> transformDistance2b(const QuipDescriptorData &main_data,
                                                                 const pugi::xml_node &distance2b_node);

        /** Builds a 3-body component (Angle3bTransformation + SquaredExpKernel) from an {@code angle_3b}
         *  descriptor. */
        static ValuePtr<GapComponent> transformAngle3b(const QuipDescriptorData &mainData,
                                                              const pugi::xml_node &angle3b_node);

        /** Builds an EAM many-body component from an {@code eam_density} descriptor, given the set of
         *  species present. */
        static ValuePtr<GapComponent> transformEam(const QuipDescriptorData &main_data,
                                                          const pugi::xml_node &eam_nodes,
                                                          const std::set<Species> &species);

        /** Selects the EAM base pair (density) function (FSGen / Polycutoff / Coscutoff / spline) from the
         *  descriptor's {@code pair_function}/{@code order}/{@code mode} fields. */
        static ValuePtr<ClusterTransformation<1, 2> > selectPairFunction(const QuipDescriptorData &main_data,
                                                                                std::optional<double> r_min,
                                                                                double prefactor);
    };
}

#endif