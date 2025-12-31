#include "ThreeBodySE.hpp"

#include <utility>

#include "utils/Utils.hpp"

namespace jgap {
    ThreeBodySE::ThreeBodySE(SpeciesTriplet idTriplet, double energyScale, Vector3 lengthScales,
                             Vector3 q, double fCut)
        : _idTriplet(std::move(idTriplet)),
          _energyScale(energyScale),
          _lengthScale(lengthScales),
          _q(q),
          _descriptorPrefactors(fCut) {
        _totalPrefactor = _descriptorPrefactors;
        _inverseThetaSq = {
            1.0 / (_lengthScale.x * _lengthScale.x),
            1.0 / (_lengthScale.y * _lengthScale.y),
            1.0 / (_lengthScale.z * _lengthScale.z)
        };
    }

    ThreeBodySE::ThreeBodySE(const DataNode &params)
        : _idTriplet(SpeciesTriplet{
              REQUIRE(params, "species_triplet")[0],
              SpeciesPair{REQUIRE(params, "species_triplet")[1], REQUIRE(params, "species_triplet")[2]}
          }),
          _energyScale(REQUIRE(params, "energy_scale")),
          _descriptorPrefactors(REQUIRE(params, "descriptor_prefactors")),
          _q({REQUIRE(params, "q")[0], REQUIRE(params, "q")[1], REQUIRE(params, "q")[2]}) {
        auto lsJson = REQUIRE(params, "length_scale");
        if (lsJson.is_number()) {
            _lengthScale = {lsJson, lsJson, lsJson};
        } else {
            _lengthScale = {lsJson[0], lsJson[1], lsJson[2]};
        }

        _totalPrefactor = _descriptorPrefactors * pow(_energyScale, 2);
        _inverseThetaSq = {
            1.0 / (_lengthScale.x * _lengthScale.x),
            1.0 / (_lengthScale.y * _lengthScale.y),
            1.0 / (_lengthScale.z * _lengthScale.z)
        };

        if (params.contains("coefficient")) {
            coefficient = params["coefficient"];
        }
    }

    DataNode ThreeBodySE::serialize() {
        DataNode res = {
            {"species_triplet", std::vector{_idTriplet.root, _idTriplet.nodes.first(), _idTriplet.nodes.second()}},
            {"length_scale", std::vector{_lengthScale.x, _lengthScale.y, _lengthScale.z}},
            {"energy_scale", _energyScale},
            {"q", std::vector{_q.x, _q.y, _q.z}},
            {"descriptor_prefactors", _descriptorPrefactors}
        };
        if (coefficient.has_value()) {
            res["coefficient"] = coefficient.value();
        }
        return res;
    }

    double ThreeBodySE::crossCovariance(const std::shared_ptr<IKernel> &other) {
        const auto otherSE = std::dynamic_pointer_cast<ThreeBodySE>(other);

        if (!otherSE) {
            JGAP_LOG_WARN("Comparing 3b SE with non 3b-SE");
            return 0.0;
        }

        return _descriptorPrefactors * otherSE->_descriptorPrefactors
               * _energyScale * otherSE->_energyScale
               * exp(-pow(_q.x - otherSE->_q.x, 2.0) / (2.0 * _lengthScale.x * otherSE->_lengthScale.x))
               * exp(-pow(_q.y - otherSE->_q.y, 2.0) / (2.0 * _lengthScale.y * otherSE->_lengthScale.y))
               * exp(-pow(_q.z - otherSE->_q.z, 2.0) / (2.0 * _lengthScale.z * otherSE->_lengthScale.z));
    }

    double ThreeBodySE::valueInternal(const Vector3 &q) const {
        return _totalPrefactor
               * exp(-pow(_q.x - q.x, 2.0) * 0.5 * _inverseThetaSq.x)
               * exp(-pow(_q.y - q.y, 2.0) * 0.5 * _inverseThetaSq.y)
               * exp(-pow(_q.z - q.z, 2.0) * 0.5 * _inverseThetaSq.z);
    }

    Vector3 ThreeBodySE::gradientInternal(const Vector3 &changingQ) const {
        return (_q - changingQ).componentMul(_inverseThetaSq) * valueInternal(changingQ);
    }
}
