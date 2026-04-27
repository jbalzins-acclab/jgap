
#include "core/descriptors/eam/EamSE.hpp"

#include <utility>

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

    Covariance EamSE::covariance(const AtomicStructure &structure,
                                 const EamKernelIndexPerSpecies &indexes) {

        double energy = 0;
        vector forces(structure.size(), Vector3{0.0, 0.0, 0.0});
        array<Vector3, 3> virials{};

        for (const auto &index : indexes) {
            energy += value(index.density);

            const double dU_drho_i = derivative(index.density);
            auto atomPosition = structure.positions[index.atAtomIndex];

            for (auto &[neighbourData, d_rho_i_dr_ij]: index.densityDerivatives) {
                const Vector3 r01 = structure.positions[neighbourData.index] + neighbourData.offset - atomPosition;
                const Vector3 f10 = r01.normalize() * d_rho_i_dr_ij * dU_drho_i;
                forces[index.atAtomIndex] += f10;
                forces[neighbourData.index] -= f10;

                // x2 since r10.x * f10.x = r01.x * f01.x // what ???
                virials[0] -= f10 * r01.x;
                virials[1] -= f10 * r01.y;
                virials[2] -= f10 * r01.z;
            }
        }

        return {energy, forces, virials};
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

    double EamSE::value(const double &density) {
        return _totalPrefactor * exp(-pow(density - _density, 2) * 0.5 * _inverseThetaSq);
    }

    double EamSE::derivative(const double &density) {
        return (_density - density) * _inverseThetaSq * value(density);
    }
}
