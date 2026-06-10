#ifndef QUIPXMLCONVERTER_HPP
#define QUIPXMLCONVERTER_HPP

#include <pugixml.hpp>

#include "core/potentials/Potential.hpp"
#include "core/potentials/gap/GapPotential.hpp"
#include "core/potentials/isolated/IsolatedAtomPotential.hpp"

namespace jgap {
    class QuipXmlConverter {
    public:
        static std::unique_ptr<Potential> transform(const pugi::xml_node& quip_potential_encoded);
        static std::unique_ptr<Potential> transformPairpot(const pugi::xml_node& quip_pairpot);
        static std::unique_ptr<Potential> transformGapParams(const pugi::xml_node& quip_gap_params);

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

        static std::unique_ptr<IsolatedAtomPotential> transformIsolatedAtomParams(
                                                        const pugi::xml_node& quip_isolated_atom_params);

        static std::unique_ptr<GapPotential> transformSparseData(const pugi::xml_node& quip_sparse_data);

        static std::unique_ptr<GapComponent> transformDistance2b(const QuipDescriptorData &main_data,
                                                                 const pugi::xml_node &distance2b_node);

        static std::unique_ptr<GapComponent> transformAngle3b(const QuipDescriptorData &mainData,
                                                              const pugi::xml_node &angle3b_node);

        static std::unique_ptr<GapComponent> transformEam(const QuipDescriptorData &main_data,
                                                          const pugi::xml_node &eam_nodes,
                                                          const std::set<Species> &species);

        static std::unique_ptr<ClusterTransformation<1, 2> > selectPairFunction(const QuipDescriptorData &main_data,
                                                                                std::optional<double> r_min,
                                                                                double prefactor);
    };
}

#endif