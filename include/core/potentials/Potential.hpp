#ifndef EXTERNALPOTENTIAL_HPP
#define EXTERNALPOTENTIAL_HPP

#include "data/BasicDataTypes.hpp"
#include "data/TabulationData.hpp"

#include <memory>
#include <nlohmann/json.hpp>

namespace jgap {
    class Potential {
    public:
        virtual ~Potential() = default;
        virtual PotentialPrediction predict(const AtomicStructure& structure) = 0;
        virtual nlohmann::json serialize() = 0;
        virtual string getType() = 0;
        virtual double getCutoff() = 0;

        virtual TabulationData tabulate(const TabulationParams& params) = 0;
    };
}

#endif
