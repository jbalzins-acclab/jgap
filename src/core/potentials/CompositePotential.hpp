#ifndef JGAP_COMPOSITEPOTENTIALFIT_HPP
#define JGAP_COMPOSITEPOTENTIALFIT_HPP

#include "Potential.hpp"
#include "../../data/atomic/AtomicStructure.hpp"
#include "core/tabulation/Tabulatable.hpp"
#include "io/Serializable.hpp"
#include "io/parse/ParserRegistry.hpp"

namespace jgap {

    class CompositePotential : public Potential, Serializable, Tabulatable {
    public:
        SETUP_PARSER_AND_SERIALIZATION(Potential, CompositePotential, composite)

        std::map<std::string, std::shared_ptr<Potential>> potentials{};

        ~CompositePotential() override = default;
        CompositePotential(const std::map<std::string, std::shared_ptr<Potential>> &potentialsMap);

        CutoffRanges getCutoff() override;

        Predictions predict(const AtomicStructure &structure) override;

        void tabulate(TabulationData& table) override;
    };
}

#endif
