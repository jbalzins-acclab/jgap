#ifndef QUIPXMLCONVERTER_HPP
#define QUIPXMLCONVERTER_HPP

#include <nlohmann/json.hpp>
#include <pugixml.hpp>

#include "core/descriptors/EamDescriptor.hpp"
#include "core/descriptors/ThreeBodyDescriptor.hpp"
#include "core/descriptors/TwoBodyDescriptor.hpp"
#include "core/potentials/IsolatedAtomPotential.hpp"
#include "core/potentials/JgapPotential.hpp"
#include "core/potentials/Potential.hpp"

using namespace std;

namespace jgap {
    class QuipXmlConverter {
    public:
        static shared_ptr<Potential> transform(pugi::xml_node quipPotential);
        static shared_ptr<Potential> transformPairpot(pugi::xml_node quipPairpot);
        static shared_ptr<Potential> transformGapParams(pugi::xml_node quipGapParams);

    private:
        struct QuipDescriptorData {
            string type;
            double delta;
            double theta;
            double cutoff;

            optional<double> rMin;
            optional<double> cutoffTransitionWidth;
            optional<string> pairFunction;
            optional<double> order;
            optional<string> mode;

            bool operator==(const QuipDescriptorData &other) const;
            bool operator<(const QuipDescriptorData &other) const;
        };

        static shared_ptr<IsolatedAtomPotential> transformIsolatedAtomParams(pugi::xml_node quipIsolatedAtomParams);

        static shared_ptr<JgapPotential> transformSparseData(pugi::xml_node quipSparseData);
        static shared_ptr<TwoBodyDescriptor> transformDistance2b(QuipDescriptorData mainData,
                                                                 vector<pugi::xml_node> distance2bNodes);
        static shared_ptr<ThreeBodyDescriptor> transformAngle3b(QuipDescriptorData mainData,
                                                                vector<pugi::xml_node> distance2bNodes);
        static shared_ptr<EamDescriptor> transformEam(QuipDescriptorData mainData,
                                                      const vector<pugi::xml_node> &eamNodes);
        static shared_ptr<EamPairFunction> selectPairFunction(QuipDescriptorData mainData,
                                                              optional<double> rMin, double prefactor);
    };
}

#endif