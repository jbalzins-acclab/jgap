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
        explicit GapPotential(const nlohmann::json& params);
        explicit GapPotential(map<string, shared_ptr<Descriptor>> descriptors)
            : _descriptors(std::move(descriptors)) {
        }
        ~GapPotential() override = default;

        PotentialPrediction predict(const AtomicStructure &structure) override;

        nlohmann::json serialize() override;
        string getType() override { return "gap"; }
        double getCutoff() override;

        TabulationData tabulate(const TabulationParams &params) override;

    private:
        map<string, shared_ptr<Descriptor>> _descriptors;
    };

    REGISTER_PARSER("gap", Potential, GapPotential);
}
#endif
