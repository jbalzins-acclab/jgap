#ifndef GAPPOTENTIAL_HPP
#define GAPPOTENTIAL_HPP

#include <string_view>
#include <utility>
#include <map>
#include <memory>
#include "data/DataNode.hpp"

#include "core/descriptors/Descriptor.hpp"
#include "core/potentials/Potential.hpp"
#include "core/tabulation/Tabulatable.hpp"
#include "io/parse/ParserRegistry.hpp"

namespace jgap {
    class GapPotential : public Potential, Serializable, Tabulatable {
    public:
        SETUP_PARSER_AND_SERIALIZATION(Potential, GapPotential, gap);

        ~GapPotential() override = default;
        GapPotential(std::map<std::string, std::shared_ptr<Descriptor>> descriptors)
            : _descriptors(std::move(descriptors)) {
        }

        Predictions predict(const AtomicStructure &structure) override;

        CutoffRanges getCutoff() override;

        void tabulate(TabulationData& table) override;

    private:
        std::map<std::string, std::shared_ptr<Descriptor>> _descriptors;
    };
}
#endif
