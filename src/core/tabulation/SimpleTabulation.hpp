#ifndef JGAP_SIMPLETABULATION_HPP
#define JGAP_SIMPLETABULATION_HPP

#include <memory>

#include "Tabulation.hpp"
#include "core/potentials/Potential.hpp"
#include "core/potentials/TabGapPotential.hpp"

namespace jgap {
    class SimpleTabulation : public Tabulation {
    public:
        SETUP_PARSER(Tabulation, SimpleTabulation, default)

        SimpleTabulation(TabulationParams defaultParams, std::optional<std::string> tableFilenamePrefix);

        std::shared_ptr<TabGapPotential> tabulate(const std::shared_ptr<Potential>& potential) override;

    protected:
        static std::pair<TabulationData, TabulationData> makeSplineTables(const std::shared_ptr<Potential> &potential,
                                                                          const TabulationParams &params);

    private:
        TabulationParams default_params_;
        std::optional<std::string> table_filename_prefix_;

        TabulationParams prepareParams(const std::shared_ptr<Potential> &potential) const;
    };
}

#endif
