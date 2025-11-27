#ifndef JGAP_SIMPLETABULATION_HPP
#define JGAP_SIMPLETABULATION_HPP

#include <memory>
#include <nlohmann/json.hpp>

#include "Tabulation.hpp"
#include "core/potentials/Potential.hpp"
#include "core/potentials/TabGapPotential.hpp"

namespace jgap {
    class SimpleTabulation : public Tabulation {
    public:
        static constexpr std::string TYPE = "default";
        static std::shared_ptr<SimpleTabulation> fromJson(const nlohmann::json& params);

        SimpleTabulation(TabulationParams defaultParams, std::optional<std::string> tableFilenamePrefix);

        std::shared_ptr<TabGapPotential> tabulate(const std::shared_ptr<Potential>& potential) override;

    protected:
        static std::pair<TabulationData, TabulationData> makeSplineTables(const std::shared_ptr<Potential> &potential,
                                                                          const TabulationParams &params);

    private:
        TabulationParams _defaultParams;
        std::optional<std::string> _tableFilenamePrefix;

        TabulationParams prepareParams(const std::shared_ptr<Potential>& potential) const;
    };
    REGISTER_PARSER(Tabulation, SimpleTabulation)
}

#endif
