#ifndef JGAPPOTENTIAL_HPP
#define JGAPPOTENTIAL_HPP

#include <string_view>
#include <utility>
#include <nlohmann/json.hpp>

#include "core/potentials/Potential.hpp"

namespace jgap {
    class JgapPotential : public Potential {
    public:
        JgapPotential(vector<shared_ptr<Descriptor>> descriptors)
            : _descriptors(descriptors) {}

        ~JgapPotential() override = default;

        PotentialPrediction predict(const AtomicStructure &structure) override {
            PotentialPrediction prediction{};
            for (const auto &descriptor : _descriptors) {
                prediction = prediction + descriptor->predict(structure);
            }
            return prediction;
        }

        nlohmann::json serialize() override {
            nlohmann::json descriptors = nlohmann::json::array();

            for (const auto &descriptor : _descriptors) {
                descriptors.push_back(descriptor->serialize());
            }

            return {
                {"potential", "jgap"},
                {"descriptors", descriptors}
            };
        }

    private:
        vector<shared_ptr<Descriptor>> _descriptors;
    };
}

#endif
