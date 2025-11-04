
#include "core/descriptors/eam/EamSE.hpp"

#include <utility>

#include "io/log/CurrentLogger.hpp"
#include "utils/Utils.hpp"

namespace jgap {

    EamSE::EamSE(Species species, double energyScale, double lengthScale, double density, optional<double> coeff)
        : _idSpecies(std::move(species)),
          _energyScale(energyScale),
          _lengthScale(lengthScale),
          _density(density) {

        _inverseThetaSq = 1.0 / (_lengthScale * _lengthScale);
        _totalPrefactor = energyScale * energyScale;

        coefficient = coeff;
    }

    EamSE::EamSE(const nlohmann::json &params)
        : _idSpecies(require(params, "species")),
          _energyScale(require(params, "energy_scale")),
          _lengthScale(require(params, "length_scale")),
          _density(require(params, "density")) {

        _totalPrefactor = _energyScale * _energyScale;
        _inverseThetaSq = 1.0 / (_lengthScale * _lengthScale);

        if (params.contains("coefficient")) {
            coefficient = params["coefficient"];
        }
    }

    nlohmann::json EamSE::serialize() {
        nlohmann::json res = {
            {"species", _idSpecies},
            {"length_scale", _lengthScale},
            {"energy_scale", _energyScale},
            {"density", _density}
        };
        if (coefficient.has_value()) {
            res["coefficient"] = coefficient.value();
        }
        return res;
    }


    double EamSE::crossCovariance(const shared_ptr<IKernel> &other) {

        const auto otherSE = std::dynamic_pointer_cast<EamSE>(other);

        if (!otherSE) {
            CurrentLogger::get()->warn("Comparing EAM SE with non EAM SE");
            return 0.0;
        }

        return _energyScale * otherSE->_energyScale
               * exp(-pow(_density - otherSE->_density, 2.0) / (2.0 * _lengthScale * otherSE->_lengthScale));
    }

    double EamSE::valueInternal(const double &density) const {
        return _totalPrefactor * exp(-pow(density - _density, 2) * 0.5 * _inverseThetaSq);
    }

    double EamSE::derivativeInternal(const double &density) const {
        return (_density - density) * _inverseThetaSq * valueInternal(density);
    }
}
