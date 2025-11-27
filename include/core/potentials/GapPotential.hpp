#ifndef GAPPOTENTIAL_HPP
#define GAPPOTENTIAL_HPP

#include <string_view>
#include <utility>
#include <nlohmann/json.hpp>

#include "core/descriptors/Descriptor.hpp"
#include "core/potentials/Potential.hpp"
#include "io/parse/ParserRegistry.hpp"

namespace jgap {
    class GapPotential : public Potential {
    public:
        static constexpr string TYPE = "gap";

        explicit GapPotential(const nlohmann::json& params);
        explicit GapPotential(std::std::map<string, shared_ptr<Descriptor>> descriptors)
            : _descriptors(std::move(descriptors)) {
        }
        ~GapPotential() override = default;

        Predictions predict(const AtomicStructure &structure) override;

        nlohmann::json serialize() override;
        string getType() override { return TYPE; }
        CutoffRanges getCutoff() override;

        void tabulate(TabulationData& table) override;

    private:
        std::map<string, shared_ptr<Descriptor>> _descriptors;
    };

    REGISTER_PARSER(Potential, GapPotential);
}
#endif
