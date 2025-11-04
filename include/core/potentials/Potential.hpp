#ifndef JGAP_POTENTIAL_HPP
#define JGAP_POTENTIAL_HPP

#include "data/TabulationData.hpp"

#include <nlohmann/json.hpp>

#include "data/AtomicStructure.hpp"
#include "data/CutoffRanges.hpp"
#include "data/PredictionData.hpp"

namespace jgap {
    class Potential {
    public:
        virtual ~Potential() = default;
        virtual Predictions predict(const AtomicStructure& structure) = 0;
        virtual nlohmann::json serialize() = 0;
        virtual string getType() = 0;
        virtual CutoffRanges getCutoff() = 0;

        virtual void tabulate(TabulationData& table) = 0;
    };
}

#endif
