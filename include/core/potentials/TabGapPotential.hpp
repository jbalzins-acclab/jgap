#ifndef JGAP_TABGAPPOTENTIAL_HPP
#define JGAP_TABGAPPOTENTIAL_HPP

#include <io/parse/ParserRegistry.hpp>

#include "Potential.hpp"
#include "core/descriptors/2b/TwoBodyDescriptor.hpp"
#include "core/descriptors/3b/ThreeBodyDescriptor.hpp"
#include "core/descriptors/eam/EamDescriptor.hpp"
#include "data/AtomicStructure.hpp"
#include "data/PredictionData.hpp"

namespace jgap {

    class TabGapPotential : public Potential {
    public:
        static constexpr string TYPE = "tabgap";

        ~TabGapPotential() override = default;
        TabGapPotential(const nlohmann::json& params);
        TabGapPotential(TabulationData splineCoefficients, const vector<string>& files);

        Predictions predict(const AtomicStructure &structure) override;

        nlohmann::json serialize() override;

        string getType() override { return TYPE; }
        CutoffRanges getCutoff() override;

        void tabulate(TabulationData& table) override;

    private:
        nlohmann::json _params;

        map<Species, double> _isolatedEnergies;
        shared_ptr<TwoBodyDescriptor> _2b;
        shared_ptr<ThreeBodyDescriptor> _3b;
        vector<shared_ptr<EamDescriptor>> _eams;

        void init(TabulationData& splineCoefficients);
    };

    REGISTER_PARSER(Potential, TabGapPotential)
}

#endif
