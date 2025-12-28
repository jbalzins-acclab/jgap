#ifndef JGAP_TABGAPPOTENTIAL_HPP
#define JGAP_TABGAPPOTENTIAL_HPP

#include <io/parse/ParserRegistry.hpp>

#include "Potential.hpp"
#include "core/descriptors/2b/TwoBodyDescriptor.hpp"
#include "core/descriptors/3b/ThreeBodyDescriptor.hpp"
#include "core/descriptors/eam/EamDescriptor.hpp"
#include "../../data/atomic/AtomicStructure.hpp"
#include "../../data/atomic/PredictionData.hpp"

namespace jgap {

    class TabGapPotential : public Potential, Serializable {
    public:
        SETUP_PARSER_AND_SERIALIZATION(Potential, TabGapPotential, tabgap)

        ~TabGapPotential() override = default;
        TabGapPotential(TabulationData spline_coefficients, const std::vector<std::string>& files);

        Predictions predict(const AtomicStructure &structure) override;

        CutoffRanges getCutoff() override;

        void tabulate(TabulationData& table) override;

    private:
        DataNode params_;

        std::map<Species, double> isolated_energies_;
        std::shared_ptr<TwoBodyDescriptor> two_body_;
        std::shared_ptr<ThreeBodyDescriptor> three_body_;
        std::vector<std::shared_ptr<EamDescriptor>> eams_;

        void init(TabulationData& spline_coefficients);
    };
}

#endif
