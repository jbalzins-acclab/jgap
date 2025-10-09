#include <utility>

#include "core/descriptors/2b/TwoBodySE.hpp"

#include "io/parse/ParserRegistry.hpp"
#include "utils/Utils.hpp"

namespace jgap {
    TwoBodySE::TwoBodySE(SpeciesPair speciesPair, double energyScale, double lengthScale, double r, double fCut,
                         optional<double> coeff)
        : _idPair(std::move(speciesPair)),
          _energyScale(energyScale),
          _lengthScale(lengthScale),
          _r(r),
          _descriptorPrefactors(fCut) {

        _totalPrefactor = _descriptorPrefactors * _energyScale * _energyScale;
        _inverseThetaSq = 1.0 / (_lengthScale * _lengthScale);

        coefficient = coeff;
    }

    TwoBodySE::TwoBodySE(const nlohmann::json &params)
        : _idPair(SpeciesPair{require(params, "species_pair")[0], require(params, "species_pair")[1]}),
          _energyScale(require(params, "energy_scale")),
          _lengthScale(require(params, "length_scale")),
          _r(require(params, "r")),
          _descriptorPrefactors(require(params, "descriptor_prefactors")) {

        _totalPrefactor = _descriptorPrefactors * _energyScale * _energyScale;
        _inverseThetaSq = 1.0 / (_lengthScale * _lengthScale);

        if (params.contains("coefficient")) {
            coefficient = params["coefficient"];
        }
    }

    nlohmann::json TwoBodySE::serialize() {
        nlohmann::json res = {
            {"species_pair", vector{_idPair.first(), _idPair.second()}},
            {"length_scale", _lengthScale},
            {"energy_scale", _energyScale},
            {"r", _r},
            {"descriptor_prefactors", _descriptorPrefactors}
        };
        if (coefficient.has_value()) {
            res["coefficient"] = coefficient.value();
        }
        return res;
    }

    Covariance TwoBodySE::covariance(const AtomicStructure &structure, const TwoBodyIndex &index) {
        double energy = 0;
        vector<Vector3> forces(structure.size(), {0, 0, 0});
        array<Vector3, 3> virials{};

        for (const auto &[atomIndex0, atomIndex1, r01, r, fCut, dCut_dr]: index) {
            // ---------------------- ENERGY --------------------------------
            auto cov = valueInternal(r) * fCut;
            if (atomIndex0 != atomIndex1)
                cov *= 2.0; // K(r_ij,)+K(r_ji)(?)
            // else *= 1.0 --- since both r_ii(+offset) and r_ii(-offset) are in index
            energy += cov;

            // ---------------------- FORCES --------------------------------
            double dE_dr = derivativeInternal(r) * fCut;

            if (fCut < 1.0) {
                dE_dr += valueInternal(r) * dCut_dr;
            }

            auto f10 = r01.normalize() * dE_dr;
            if (atomIndex0 != atomIndex1) {
                f10 *= 2.0;
                forces[atomIndex0] += f10;
                forces[atomIndex1] -= f10;
            }

            virials[0] -= f10 * r01.x;
            virials[1] -= f10 * r01.y;
            virials[2] -= f10 * r01.z;
        }

        return {energy, forces, virials};
    }

    double TwoBodySE::value(const TwoBodyDescriptorData &r) {
        return valueInternal(r.r) * r.fCut;
    }

    double TwoBodySE::crossCovariance(const shared_ptr<IKernel> &other) {
        const auto otherSE = std::dynamic_pointer_cast<TwoBodySE>(other);

        if (!otherSE) {
            CurrentLogger::get()->warn("Comparing 2b SE with non 2b-SE");
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
