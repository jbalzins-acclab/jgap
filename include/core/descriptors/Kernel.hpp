#ifndef JGAP_KERNEL_HPP
#define JGAP_KERNEL_HPP

#include "data/Vector3.hpp"
#include <optional>
#include <string>
#include <memory>
#include <array>
#include <vector>

#include "../../data/atomic/AtomicStructure.hpp"
#include "data/DataNode.hpp"
#include "../../data/atomic/PredictionData.hpp"
#include "io/Serializable.hpp"
#include "io/log/CurrentLogger.hpp"

namespace jgap {
    class IKernel {
    public:
        std::optional<double> coefficient{};

        virtual ~IKernel() = default;
        virtual double crossCovariance(const std::shared_ptr<IKernel> &kernel) = 0;
    };

    template<class TFilter, class TIndex, class TDescriptorData>
    class Kernel : public IKernel {
    public:
        ~Kernel() override = default;

        virtual Covariance covariance(const AtomicStructure&, const TIndex&) = 0;
        virtual double value(const TDescriptorData&) = 0;
        virtual TFilter getFilter() = 0;

        Predictions predict(const AtomicStructure& structure, const TIndex& indexes);
    };

    template<class TFilter, class TIndex, class TDescriptorData>
    Predictions Kernel<TFilter, TIndex, TDescriptorData>::predict(const AtomicStructure& structure,
                                                                  const TIndex& indexes) {
        if (!coefficient.has_value()) {
            JGAP_LOG_AND_THROW("Kernel coefficient not set");
        }

        const Covariance structure_covariance = covariance(structure, indexes);

        Predictions result{};
        result.energy = *coefficient * structure_covariance.total;
        result.virials = std::array{
            structure_covariance.virials[0] * (*coefficient),
            structure_covariance.virials[1] * (*coefficient),
            structure_covariance.virials[2] * (*coefficient)
        };
        result.forces = std::vector<Vector3>();
        result.forces->reserve(structure.size());
        for (const auto &force: structure_covariance.forces) {
            result.forces->push_back(force * (*coefficient));
        }

        return result;
    }
}

#endif
