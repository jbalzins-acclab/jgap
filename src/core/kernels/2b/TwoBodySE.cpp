#include <utility>
#include <memory>
#include <optional>

#include "core/descriptors/2b/TwoBodySE.hpp"

#include "io/log/CurrentLogger.hpp"
#include "io/parse/ParserRegistry.hpp"
#include "utils/Utils.hpp"

namespace jgap {
    TwoBodySE::TwoBodySE(SpeciesPair speciesPair, double energyScale, double lengthScale, double r, double fCut,
                         std::optional<double> coeff)
        : _idPair(std::move(speciesPair)),
          _energyScale(energyScale),
          _lengthScale(lengthScale),
          _r(r),
          _descriptorPrefactors(fCut) {

        _totalPrefactor = _descriptorPrefactors * _energyScale * _energyScale;
        _inverseThetaSq = 1.0 / (_lengthScale * _lengthScale);

        coefficient = coeff;
    }

    TwoBodySE::TwoBodySE(const jgap::DataNode &params)
        : _idPair(SpeciesPair{
              static_cast<std::string>(REQUIRE(params, "species_pair")[0]),
              static_cast<std::string>(REQUIRE(params, "species_pair")[1])
          }),
          _energyScale(REQUIRE(params, "energy_scale").asDouble()),
          _lengthScale(REQUIRE(params, "length_scale").asDouble()),
          _r(REQUIRE(params, "r").asDouble()),
          _descriptorPrefactors(REQUIRE(params, "descriptor_prefactors").asDouble()) {

        _totalPrefactor = _descriptorPrefactors * _energyScale * _energyScale;
        _inverseThetaSq = 1.0 / (_lengthScale * _lengthScale);

        if (params.contains("coefficient")) {
            coefficient = params.value("coefficient", 0.0);
        }
    }

    jgap::DataNode TwoBodySE::serialize() {
        jgap::DataNode res = {
            {"species_pair", std::vector{_idPair.first(), _idPair.second()}},
            {"length_scale", _lengthScale},
            {"energy_scale", _energyScale},
            {"r", _r},
            {"descriptor_prefactors", _descriptorPrefactors}
        };
        if (coefficient.has_value()) {
            res["coefficient"] = jgap::DataNode(coefficient.value());
        }
        return res;
    }

    double TwoBodySE::crossCovariance(const std::shared_ptr<IKernel> &other) {
        const auto otherSE = std::dynamic_pointer_cast<TwoBodySE>(other);

        if (!otherSE) {
            JGAP_LOG_WARN("Comparing 2b-SE with non 2b-SE");
            return 0.0;
        }

        return _descriptorPrefactors * otherSE->_descriptorPrefactors
               * _energyScale * otherSE->_energyScale
               * exp(-pow(_r - otherSE->_r, 2.0) / (2.0 * _lengthScale * otherSE->_lengthScale));
    }

    double TwoBodySE::valueInternal(const double &r) const {
        return _totalPrefactor * exp(-pow(r - _r, 2.0) * 0.5 * _inverseThetaSq);
    }

    double TwoBodySE::derivativeInternal(const double &changingR) const {
        return (_r - changingR) * _inverseThetaSq * valueInternal(changingR);
    }
}
