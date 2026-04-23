#ifndef ISOLATEDATOMPOTENTIAL_HPP
#define ISOLATEDATOMPOTENTIAL_HPP

#include <map>
#include "Potential.hpp"
#include "core/tabulation/Tabulatable.hpp"
#include "io/Serializable.hpp"
#include "io/parse/ParserRegistry.hpp"

namespace jgap {
    class IsolatedAtomPotential : public Potential, Serializable, Tabulatable {
    public:
        SETUP_PARSER_AND_SERIALIZATION(Potential, IsolatedAtomPotential, isolated_atom);

        ~IsolatedAtomPotential() override = default;
        IsolatedAtomPotential(const std::map<Species, double>& isolated_atom_energies, bool error_on_unknown);

        CutoffRanges getCutoff() override { return {}; }

        Predictions predict(const AtomicStructure& structure) override;

        void tabulate(TabulationData& table) override;

    private:
        bool error_on_unknown_species_;
        std::map<Species, double> isolated_energies_;
    };
}

#endif
