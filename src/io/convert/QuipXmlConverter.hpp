#ifndef JGAP_QUIPXMLCONVERTER_HPP
#define JGAP_QUIPXMLCONVERTER_HPP

#include <pugixml.hpp>

#include <filesystem>
#include <string>

#include "core/potentials/Potential.hpp"
#include "core/potentials/gap/GapPotential.hpp"
#include "core/potentials/isolated/IsolatedAtomPotential.hpp"
#include "core/transform/eam/EamPairFunction.hpp"

namespace jgap {
    class QuipXmlConverter {
    public:
        static ValuePtr<Potential> transform(const std::filesystem::path& xml_filename);

    private:
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

        static ValuePtr<Potential> transform(const pugi::xml_node& quip_potential_encoded,
                                              const std::filesystem::path& base_dir);

        static ValuePtr<Potential> transformPairpot(const pugi::xml_node& quip_pairpot);

        static ValuePtr<Potential> transformGapParams(const pugi::xml_node& quip_gap_params,
                                                      const std::filesystem::path& base_dir);


        /** Reads the {@code GAP_data} e0 entries into an {@link IsolatedAtomPotential}. */
        static IsolatedAtomPotential transformIsolatedAtomParams(
                                                        const pugi::xml_node& quip_isolated_atom_params);

        /** Builds the GAP components from the {@code gpSparse} block (one per {@code gpCoordinates}).
         *  @param base_dir directory the referenced {@code sparseX} filenames are resolved against. */
        static GapPotential transformSparseData(const pugi::xml_node& quip_sparse_data,
                                                const std::filesystem::path& base_dir);

        /**
         * Builds a 2-body component (TwoBodyTransformation + SquaredExpKernel) from a {@code distance_2b}
         * descriptor, reading sparse points from the referenced {@code sparseX} file (resolved against
         * {@code base_dir}) and the per-point {@code alpha} coefficients.
         */
        static ValuePtr<GapComponent> transformDistance2b(const QuipDescriptorData &main_data,
                                                                 const pugi::xml_node &distance2b_node,
                                                                 const std::filesystem::path& base_dir);

        /** Builds a 3-body component (Angle3bTransformation + SquaredExpKernel) from an {@code angle_3b}
         *  descriptor. */
        static ValuePtr<GapComponent> transformAngle3b(const QuipDescriptorData &mainData,
                                                              const pugi::xml_node &angle3b_node,
                                                              const std::filesystem::path& base_dir);

        /** Builds an EAM many-body component from an {@code eam_density} descriptor, given the set of
         *  species present. */
        static ValuePtr<GapComponent> transformEam(const QuipDescriptorData &main_data,
                                                          const pugi::xml_node &eam_nodes,
                                                          const std::set<Species> &species,
                                                          const std::filesystem::path& base_dir);

        /** Resolves a {@code sparseX_filename} against {@code base_dir}: absolute paths and the
         *  empty-{@code base_dir} case are returned unchanged, otherwise joined onto {@code base_dir}. */
        static std::string resolveSparseX(const std::filesystem::path& base_dir, const std::string& filename);

        /** Selects the EAM base pair (density) function (FSGen / Polycutoff / Coscutoff / spline) from the
         *  descriptor's {@code pair_function}/{@code order}/{@code mode} fields. */
        static ValuePtr<EamPairFunction> selectPairFunction(const QuipDescriptorData &main_data,
                                                                                std::optional<double> r_min,
                                                                                double prefactor);
    };
}

#endif