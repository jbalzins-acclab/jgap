#ifndef QUIPXMLCONVERTER_HPP
#define QUIPXMLCONVERTER_HPP

#include <pugixml.hpp>

#include "core/descriptors/eam/EamDescriptor.hpp"
#include "core/descriptors/3b/ThreeBodyDescriptorFinder.hpp"
#include "core/descriptors/2b/TwoBodyDescriptor.hpp"
#include "core/potentials/IsolatedAtomPotential.hpp"
#include "core/potentials/GapPotential.hpp"
#include "core/potentials/Potential.hpp"

namespace jgap {
    class QuipXmlConverter {
    public:
        static std::shared_ptr<Potential> transform(pugi::xml_node quipPotential);
        static std::shared_ptr<Potential> transformPairpot(pugi::xml_node quipPairpot);
        static std::shared_ptr<Potential> transformGapParams(pugi::xml_node quipGapParams);

    private:
        struct QuipDescriptorData {
            std::string type;
            double delta;
            double theta;
            double cutoff;

            std::optional<double> rMin;
            std::optional<double> cutoffTransitionWidth;
            std::optional<std::string> pairFunction;
            std::optional<double> order;
            std::optional<std::string> mode;

            bool operator==(const QuipDescriptorData &other) const;
            bool operator<(const QuipDescriptorData &other) const;
        };

        static std::shared_ptr<IsolatedAtomPotential> transformIsolatedAtomParams(
                                                        pugi::xml_node quipIsolatedAtomParams);

        static std::shared_ptr<GapPotential> transformSparseData(pugi::xml_node quipSparseData);
        static std::shared_ptr<TwoBodyDescriptor> transformDistance2b(QuipDescriptorData mainData,
                                                                 std::vector<pugi::xml_node> distance2bNodes);
        static std::shared_ptr<ThreeBodyDescriptorFinder> transformAngle3b(QuipDescriptorData mainData,
                                                                std::vector<pugi::xml_node> distance2bNodes);
        static std::shared_ptr<EamDescriptor> transformEam(QuipDescriptorData mainData,
                                                      const std::vector<pugi::xml_node> &eamNodes);
        static std::shared_ptr<EamPairFunction> selectPairFunction(QuipDescriptorData mainData,
                                                              std::optional<double> rMin, double prefactor);
    };
}

#endif