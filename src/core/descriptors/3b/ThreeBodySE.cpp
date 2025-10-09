#include <utility>

#include "core/descriptors/3b/ThreeBodySE.hpp"

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

    ThreeBodySE::ThreeBodySE(const nlohmann::json &params)
        : _idTriplet(SpeciesTriplet{
            require(params, "species_triplet")[0],
            SpeciesPair{require(params, "species_triplet")[1], require(params, "species_triplet")[2]}
          }),
          _energyScale(require(params, "energy_scale")),
          _descriptorPrefactors(require(params, "descriptor_prefactors")),
          _q({require(params, "q")[0], require(params, "q")[1], require(params, "q")[2]}) {

        auto lsJson = require(params, "length_scale");
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

    nlohmann::json ThreeBodySE::serialize() {
        nlohmann::json res = {
            {"species_triplet", vector{_idTriplet.root, _idTriplet.nodes.first(), _idTriplet.nodes.second()}},
            {"length_scale", vector{_lengthScale.x, _lengthScale.y, _lengthScale.z}},
            {"energy_scale", _energyScale},
            {"q", vector{_q.x, _q.y, _q.z}},
            {"descriptor_prefactors", _descriptorPrefactors}
        };
        if (coefficient.has_value()) {
            res["coefficient"] = coefficient.value();
        }
        return res;
    }

    double ThreeBodySE::crossCovariance(const shared_ptr<IKernel> &other) {
        const auto otherSE = std::dynamic_pointer_cast<ThreeBodySE>(other);

        if (!otherSE) {
            CurrentLogger::get()->warn("Comparing 3b SE with non 3b-SE");
            return 0.0;
        }

        return _descriptorPrefactors * otherSE->_descriptorPrefactors
               * _energyScale * otherSE->_energyScale
               * exp(-pow(_q.x - otherSE->_q.x, 2.0) / (2.0 * _lengthScale.x * otherSE->_lengthScale.x))
               * exp(-pow(_q.y - otherSE->_q.y, 2.0) / (2.0 * _lengthScale.y * otherSE->_lengthScale.y))
               * exp(-pow(_q.z - otherSE->_q.z, 2.0) / (2.0 * _lengthScale.z * otherSE->_lengthScale.z));
    }

    Covariance ThreeBodySE::covariance(const AtomicStructure &structure,
                                       const ThreeBodyIndex &index) {
        double energy = 0;
        vector forces(structure.size(), Vector3{0, 0, 0});
        array<Vector3, 3> virials{};

        for (auto [atomIndexes,
                   r_ij, grad_rij_wrt_rj,
                   fCut01, fCut02, dfCut_dr_01, dfCut_dr_02,
                   q, dq_k_dr_ij]: index) {
            // ---------------------- ENERGY --------------------------------
            energy += valueInternal(q)
                        * fCut01 * fCut02
                        * 2.0/*q_ijk + q_ikj*/;

            // ---------------------- FORCES --------------------------------
            // e.g. dU/dq
            Vector3 gradUWrtQ = gradientInternal(q);

            // dU/d{r01, r02, r12} = dU/dq * dq/dr_ij
            Vector3 gradUWrtDistances = Vector3{
                gradUWrtQ.dot(dq_k_dr_ij[0]),
                gradUWrtQ.dot(dq_k_dr_ij[1]),
                gradUWrtQ.dot(dq_k_dr_ij[2])
            } * fCut01 * fCut02;

            // Product rule with cutoffs
            if (fCut01 < 1.0) {
                gradUWrtDistances.x += valueInternal(q) * dfCut_dr_01 * fCut02;
            }
            if (fCut02 < 1.0) {
                gradUWrtDistances.y += valueInternal(q) * fCut01 * dfCut_dr_02;
            }
            gradUWrtDistances *= 2.0/*q_ijk + q_ikj */;

            // chain rule x2 ( remember: gradWrtDistances = d/d{r01, r02, r12} )

            const Vector3 f10 = grad_rij_wrt_rj[0] * gradUWrtDistances.x;
            forces[atomIndexes[0]] += f10;
            forces[atomIndexes[1]] -= f10;
            virials[0] -= f10 * r_ij[0].x;
            virials[1] -= f10 * r_ij[0].y;
            virials[2] -= f10 * r_ij[0].z;

            const Vector3 f20 = grad_rij_wrt_rj[1] * gradUWrtDistances.y;
            forces[atomIndexes[0]] += f20;
            forces[atomIndexes[2]] -= f20;
            virials[0] -= f20 * r_ij[1].x;
            virials[1] -= f20 * r_ij[1].y;
            virials[2] -= f20 * r_ij[1].z;

            const Vector3 f21 = grad_rij_wrt_rj[2] * gradUWrtDistances.z;
            forces[atomIndexes[1]] += f21;
            forces[atomIndexes[2]] -= f21;
            virials[0] -= f21 * r_ij[2].x;
            virials[1] -= f21 * r_ij[2].y;
            virials[2] -= f21 * r_ij[2].z;
        }

        return {energy, forces, virials};
    }

    double ThreeBodySE::value(const ThreeBodyDescriptorData &t) {
        return valueInternal(t.q) * t.fCut;
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
