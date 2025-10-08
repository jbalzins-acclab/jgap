#ifndef KERNEL_HPP
#define KERNEL_HPP

#include "data/BasicDataTypes.hpp"
#include <nlohmann/json.hpp>

namespace jgap {

    class IKernel {
    public:
        optional<double> coefficient{};

        virtual ~IKernel() = default;
        virtual string getType() = 0;
        virtual nlohmann::json serialize() = 0;

        virtual double crossCovariance(const shared_ptr<IKernel>& kernel) = 0;
    };

    template<class TFilter, class TIndex, class TDescriptorData>
    class Kernel : public IKernel {
    public:

        ~Kernel() override = default;

        virtual Covariance covariance(const AtomicStructure&, const TIndex&) = 0;
        virtual double value(const TDescriptorData&) = 0;
        virtual TFilter getFilter() = 0;

        PotentialPrediction predict(const AtomicStructure& structure, const TIndex &indexes);
    };

    template<class TFilter, class TIndex, class TDescriptorData>
    PotentialPrediction Kernel<TFilter, TIndex, TDescriptorData>::predict(const AtomicStructure &structure,
                                                                          const TIndex &indexes) {
        if (!coefficient.has_value()) {
            CurrentLogger::get()->logAndThrow("Coefficient not set in: {}", serialize().dump());
        }

        const Covariance structureCovariance = covariance(structure, indexes);

        PotentialPrediction result{};
        result.energy = *coefficient * structureCovariance.total;
        result.virials = array{
            structureCovariance.virials[0] * (*coefficient),
            structureCovariance.virials[1] * (*coefficient),
            structureCovariance.virials[2] * (*coefficient)
        };
        result.forces = vector<Vector3>();
        result.forces->reserve(structure.size());
        for (const auto &force : structureCovariance.forces) {
            result.forces->push_back(force * (*coefficient));
        }

        return result;
    }
}

#endif
