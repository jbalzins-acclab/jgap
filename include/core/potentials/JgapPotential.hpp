#ifndef JGAPPOTENTIAL_HPP
#define JGAPPOTENTIAL_HPP

#include <string_view>
#include <utility>
#include <nlohmann/json.hpp>

#include "core/descriptors/Descriptor.hpp"
#include "core/potentials/Potential.hpp"
#include "io/parse/ParserRegistry.hpp"

namespace jgap {
    class JgapPotential : public Potential {
    public:
        explicit JgapPotential(const nlohmann::json& params);
        explicit JgapPotential(map<string, shared_ptr<Descriptor>> descriptors)
            : _descriptors(std::move(descriptors)) {
        }
        ~JgapPotential() override = default;

        PotentialPrediction predict(const AtomicStructure &structure) override;

        nlohmann::json serialize() override;
        string getType() override { return "jgap"; }
        double getCutoff() override;

        TabulationData tabulate(const TabulationParams &params) override;

    private:
        map<string, shared_ptr<Descriptor>> _descriptors;
    };

    REGISTER_PARSER("jgap", Potential, JgapPotential);
}
#endif
