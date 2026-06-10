#ifdef JGAP_TABGAPPOTENTIAL_HPP
#define JGAP_TABGAPPOTENTIAL_HPP

#include <io/parse/ParserRegistry.hpp>

#include "core/potentials/Potential.hpp"

namespace jgap {

    class TabGapPotential : public Potential, Serializable {
    public:

        TabGapPotential(TabulationData spline_coefficients, const std::vector<std::string>& files);

        Predictions predict(const AtomicStructure &structure) override;

        CutoffRanges getCutoff() override;

    private:
        DataNode params_;

        std::map<Species, double> isolated_energies_;
        std::shared_ptr<TwoBodyDescriptor> two_body_;
        std::shared_ptr<ThreeBodyDescriptorFinder> three_body_;
        std::vector<std::shared_ptr<EamDescriptor>> eams_;

        void init(TabulationData& spline_coefficients);
    };
}

#endif
