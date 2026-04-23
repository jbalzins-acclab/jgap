
#include "core/descriptors/eam/EamSE.hpp"

#include <utility>

#include "io/log/CurrentLogger.hpp"
#include "utils/Utils.hpp"

namespace jgap {

    EamSE::EamSE(Species species, double energyScale, double lengthScale, double density, std::optional<double> coeff)
        : id_species_(std::move(species)),
          energy_scale_(energyScale),
          length_scale_(lengthScale),
          density_(density) {

        inverse_theta_sq = 1.0 / (length_scale_ * length_scale_);
        total_prefactor_ = energyScale * energyScale;

        coefficient = coeff;
    }

    EamSE::EamSE(const DataNode &params)
        : id_species_(require(params, "species")),
          energy_scale_(require(params, "energy_scale")),
          length_scale_(require(params, "length_scale")),
          density_(require(params, "density")) {

        total_prefactor_ = energy_scale_ * energy_scale_;
        inverse_theta_sq = 1.0 / (length_scale_ * length_scale_);

        if (params.contains("coefficient")) {
            coefficient = params["coefficient"];
        }
    }

    DataNode EamSE::serialize() {
        DataNode res = {
            {"species", id_species_},
            {"length_scale", length_scale_},
            {"energy_scale", energy_scale_},
            {"density", density_}
        };
        if (coefficient.has_value()) {
            res["coefficient"] = coefficient.value();
        }
        return res;
    }


    double EamSE::crossCovariance(const std::shared_ptr<IKernel> &other) {

        const auto otherSE = std::dynamic_pointer_cast<EamSE>(other);

        if (!otherSE) {
            JGAP_LOG_WARN("Comparing EAM SE with non EAM SE");
            return 0.0;
        }

        return energy_scale_ * otherSE->energy_scale_
               * exp(-pow(density_ - otherSE->density_, 2.0) / (2.0 * length_scale_ * otherSE->length_scale_));
    }

    double EamSE::valueInternal(const double &density) const {
        return total_prefactor_ * exp(-pow(density - density_, 2) * 0.5 * inverse_theta_sq);
    }

    double EamSE::derivativeInternal(const double &density) const {
        return (density_ - density) * inverse_theta_sq * valueInternal(density);
    }
}
