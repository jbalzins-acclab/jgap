#ifdef JGAP_TABULATION_HPP
#define JGAP_TABULATION_HPP

#include <memory>

#include "core/potentials/Potential.hpp"
#include "core/potentials/tabgap/TabGapPotential.hpp"

namespace jgap {
    class Tabulation  {
    public:

        Tabulation(TabulationParams defaultParams, std::optional<std::string> tableFilenamePrefix);

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
