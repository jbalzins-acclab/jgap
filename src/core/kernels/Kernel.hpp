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
#include "core/descriptors/Descriptor.hpp"
#include "io/Serializable.hpp"
#include "io/log/CurrentLogger.hpp"
#include "io/parse/ParserRegistry.hpp"

namespace jgap {
    class IKernel {
    public:
        std::optional<double> coefficient{};

        virtual ~IKernel() = default;
        virtual double crossCovariance(const std::shared_ptr<IKernel> &kernel) = 0;
    };

    template<size_t N_DIMENSIONS, size_t N_ATOMS>
    class RealKernel : public IKernel {
    public:
        std::array<double, N_DIMENSIONS> sparse_point;
        double sparse_f_cut;

        RealKernel(std::array<double, N_DIMENSIONS> sparse_point, double f_cut)
            : sparse_point(sparse_point), sparse_f_cut(f_cut) {
        }

        Predictions covariance(
            const AtomicStructure& structure,
            const std::vector<Descriptor<N_DIMENSIONS, N_ATOMS>>& descriptors
            ) {

            double energy = 0;
            auto forces = std::vector(structure.size(), Vector3());
            Virials virials{};

            for (const auto &descriptor: descriptors) {
                const double val = value(descriptor);
                energy += val * descriptor.f_cut;

                const auto gradU_wrt_q = gradient(descriptor); // g[i] = dU/dq_i

                for (size_t dim = 0; dim < N_DIMENSIONS; dim++) {
                    virials += descriptor.virials[dim] * gradU_wrt_q[dim];
                }

                for (const auto& gradQ_wrt_ri: descriptor.gradients) {
                    for (size_t dim = 0; dim < N_DIMENSIONS; dim++) {
                        forces[gradQ_wrt_ri[dim].wrt_atom_index] -= gradQ_wrt_ri[dim].value * gradU_wrt_q[dim];
                    }
                }
            }

            return { energy, virials, forces };
        }

        virtual double value(std::array<double, N_DIMENSIONS> q) = 0;
        virtual std::array<double, N_DIMENSIONS> gradient(std::array<double, N_DIMENSIONS> wrt_q) = 0;
    };

    template<size_t N_DIMENSIONS, size_t N_GRADIENTS>
    class SquaredExpKernel : public RealKernel<N_DIMENSIONS, N_GRADIENTS> {
    public:
        using RealKernel<N_DIMENSIONS, N_GRADIENTS>::sparse_point;

        //static constexpr std::string TYPE_ID = "squared_exp";
        //static std::shared_ptr<RealKernel<N_DIMENSIONS, N_GRADIENTS>> fromDataNode(const DataNode &params);

        SquaredExpKernel(const double energy_scale,
                         const std::array<double, N_DIMENSIONS>& length_scales,
                         const std::array<double, N_DIMENSIONS>& sparse_point,
                         const double sparse_f_cut)
            : RealKernel<N_DIMENSIONS, N_GRADIENTS>(sparse_point, sparse_f_cut),
                energy_scale_(energy_scale), length_scales_(length_scales) {

            prefactor_ = energy_scale_ * energy_scale_;
            for (size_t dim = 0; dim < N_DIMENSIONS; dim++) {
                length_scales_squared_[dim] = length_scales_[dim] * length_scales_[dim];
            }
        }

        double crossCovariance(const std::shared_ptr<IKernel> &kernel) override {
            return 0;
        }

        double value(std::array<double, N_DIMENSIONS> q) override {
            double exp_argument = 0;
            for (size_t dim = 0; dim < N_DIMENSIONS; dim++) {
                exp_argument += length_scales_squared_[dim] * (q[dim] - sparse_point[dim]);
            }
            return prefactor_ * exp(-0.5 * exp_argument);
        }

        std::array<double, N_DIMENSIONS> gradient(std::array<double, N_DIMENSIONS> wrt_q) override {

            double val = value(wrt_q);

            std::array<double, N_DIMENSIONS> result{};
            for (size_t dim = 0; dim < N_DIMENSIONS; dim++) {
                result[dim] = val * (sparse_point[dim] - wrt_q[dim]) / length_scales_squared_[dim];
            }

            return result;
        }

    private:
        double energy_scale_;
        std::array<double, N_DIMENSIONS> length_scales_;

        double prefactor_;
        std::array<double, N_DIMENSIONS> length_scales_squared_;
    };

    /*
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
    }*/
}

#endif
