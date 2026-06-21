#ifndef JGAP_TABGAPIO_HPP
#define JGAP_TABGAPIO_HPP

#include <iosfwd>
#include <memory>
#include <optional>
#include <vector>
#include <highfive/H5File.hpp>
#include <highfive/H5Group.hpp>

#include "core/potentials/tabgap/TabGapPotential.hpp"
#include "core/potentials/tabgap/components/TabGapComponent.hpp"
#include "core/potentials/tabgap/components/ThreeBodyTGComponent.hpp"
#include "core/potentials/tabgap/components/TwoBodyTGComponent.hpp"
#include "../../core/ValuePtr.hpp"

namespace jgap {

    using Filenames = std::vector<std::string>;

    /**
     * Reads and writes the externally-defined tabGAP potential format.
     *
     * Layout of a tabGAP {@code .h5}:
     * <ul>
     *   <li>{@code comment1}/{@code comment2}: free-form metadata ("UNITS: metal", "pair_style tabgap").</li>
     *   <li>{@code e0}: group of isolated-atom energies, one attribute per element symbol plus "Nelements".</li>
     *   <li>{@code npots}: [n_2b, n_3b] pair/triplet counts (n_2b is reported as 0 when the pair part is
     *       folded into the EAM .eam.fs files instead).</li>
     *   <li>{@code <i>-<j>}: one group per 2-body pair (see write2b for the grid encoding).</li>
     *   <li>{@code <i>-<j>-<k>}: one group per 3-body triplet (see write3b).</li>
     *   <li>{@code eam_files}: group of embedded LAMMPS .eam.fs file contents (one string dataset each),
     *       present only when the potential has an EAM part.</li>
     *   <li>{@code jgap}: bool attribute, true when written by jgap. Marks the file as carrying the
     *       embedded "eam_files" so they can be used as a fallback.</li>
     * </ul>
     *
     * The EAM part (embedding + density functions, and the pair potentials when EAM is present) lives in
     * LAMMPS .eam.fs files. {@link TabGapIO#write} emits them as separate files next to the .h5 AND embeds their
     * contents under "eam_files"; {@link TabGapIO#read} prefers separately-supplied .eam.fs files, but when none
     * are given and the .h5 is jgap=true it falls back to the embedded copies.
     *
     * The HDF5 core ({@link TabGapIO#writeToGroup} / {@link TabGapIO#readFromGroup}) works on a plain HighFive group, so it
     * can target either a standalone file or a node inside a larger serialization (see
     * TabGapPotentialSerialization, which reuses it and never emits separate .eam.fs files).
     */
    class TabGapIO {
    public:
        /**
         * Reads a tabGAP potential from the given files (one .h5 plus any number of .eam.fs files).
         *
         * @param filenames the files to read; exactly one must be the .h5. When no .eam.fs files are
         *                   supplied and the .h5 is jgap=true, the EAM part is taken from the embedded
         *                   "eam_files" instead.
         * @return the assembled potential.
         */
        static TabGapPotential read(const Filenames& filenames);

        /**
         * Writes a tabGAP potential to disk: an "<prefix>.tabgap.h5" plus, for the EAM part, one or more
         * "<prefix>[#i].eam.fs" files. The .eam.fs contents are also embedded inside the .h5.
         *
         * @param potential the potential to write.
         * @param output_filename_prefix prefix for the produced filenames.
         * @return the filenames written (the .h5 first, then any .eam.fs files).
         */
        static Filenames write(const TabGapPotential& potential, const std::string &output_filename_prefix);

        /**
         * Writes the whole potential into {@code root} (e0, 2b/3b groups, embedded EAM .eam.fs contents
         * under an "eam_files" group, and a jgap=true marker). Writes NO separate .eam.fs files.
         *
         * @param root the HDF5 group to write into.
         * @param potential the potential to write.
         * @return the generated .eam.fs file contents, so the file-based writer can also dump them to disk
         *         (empty when the potential has no EAM part).
         */
        static std::vector<std::string> writeToGroup(HighFive::Group& root, const TabGapPotential& potential);

        /**
         * Reads e0 and the 2b/3b groups from {@code root} into {@code pot}.
         *
         * @param root the HDF5 group to read from.
         * @param pot the potential to append the read components to.
         * @param read_embedded_eam_fs when true, the EAM part is taken from the embedded "eam_files"
         *                              datasets rather than from separate .eam.fs files.
         */
        static void readFromGroup(const HighFive::Group& root, TabGapPotential& pot, bool read_embedded_eam_fs);

        /**
         * Builds a complete TabGapPotential from a single group (with counts recomputed). Convenience for
         * callers that have just one group, e.g. TabGapPotentialSerialization.
         *
         * @param root the HDF5 group to read from.
         * @param read_embedded_eam_fs see {@link TabGapIO#readFromGroup}.
         * @return the assembled potential.
         */
        static TabGapPotential fromGroup(const HighFive::Group& root, bool read_embedded_eam_fs);

    private:
        static void write2b(HighFive::Group& root, const TwoBodyTGComponent& component);
        static void write3b(HighFive::Group& root, const ThreeBodyTGComponent& component);

        static std::string useSomeComponentsAndGenerateEamFs(const std::vector<Species>& all_species,
                                      std::map<SpeciesSet<2, FullSymmetry>, const TwoBodyTGComponent*>& pair_pots,
                                      std::multimap<Species, const EamTGComponent*>& eam_components);

        static void parseEamFs(std::istream& in, TabGapPotential& pot);
    };
}

#endif