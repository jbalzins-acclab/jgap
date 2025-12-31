#ifndef JGAP_POTENTIAL_HPP
#define JGAP_POTENTIAL_HPP

#include "../../data/tabulation/TabulationData.hpp"

#include <string>
#include <functional>

#include "../../data/atomic/AtomicStructure.hpp"
#include "data/CutoffRanges.hpp"
#include "../../data/atomic/PredictionData.hpp"
#include "data/DataNode.hpp"

namespace jgap {
    class Potential {
    public:
        virtual ~Potential() = default;

        virtual Predictions predict(const AtomicStructure& structure) = 0;
        virtual CutoffRanges getCutoff() = 0;

        virtual void tabulate(TabulationData& table) = 0;
    };
}

#endif
