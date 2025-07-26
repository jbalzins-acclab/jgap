#include "core/descriptors/kernels/TwoBodySE.hpp"

#include "io/parse/ParserRegistry.hpp"
#include "utils/Utils.hpp"

namespace jgap {
    TwoBodySE::TwoBodySE(const nlohmann::json &params) {
        _lengthScale = params["length_scale"].get<double>();
        _energyScaleSquared = pow(params["energy_scale"].get<double>(), 2);
        _inverse2ThetaSq = 1.0 / (2.0 * _lengthScale * _lengthScale);
        _inverseThetaSq = 1.0 / (_lengthScale * _lengthScale);
    }

    nlohmann::json TwoBodySE::serialize() {
        return {
            {"length_scale", _lengthScale},
            {"energy_scale", sqrt(_energyScaleSquared)}
        };
    }

    Covariance TwoBodySE::covariance(const AtomicStructure &structure,
                                     const TwoBodyKernelIndex &indexes,
                                     const TwoBodyDescriptorData &rSparse) {
        double energy = 0;
        vector<Vector3> forces(structure.size(), {0, 0, 0});
        array<Vector3, 3> virials{};

        for (const TwoBodyKernelIndexEntity &index: indexes) {
            // ---------------------- ENERGY --------------------------------
            auto cov = covarianceNoCutoffs(rSparse.r, index.r) * rSparse.fCut * index.fCut;
            if (index.atomIndex0 != index.atomIndex1) cov *= 2.0; // K(r_ij,)+K(r_ji)(?)
            energy += cov;

            // ---------------------- FORCES --------------------------------
            double dE_dr = derivativeNoCutoffs(index.r, rSparse.r) * index.fCut * rSparse.fCut;

            if (index.fCut < 1.0) {
                dE_dr += covarianceNoCutoffs(rSparse.r, index.r) * index.dCut_dr * rSparse.fCut;
            }

            auto f10 = index.r01.normalize() * dE_dr;
            if (index.atomIndex0 != index.atomIndex1) {
                f10 *= 2.0;
                forces[index.atomIndex0] += f10;
                forces[index.atomIndex1] -= f10;
            }

            virials[0] -= f10 * index.r01.x;
            virials[1] -= f10 * index.r01.y;
            virials[2] -= f10 * index.r01.z;
        }

        return {energy, forces, virials};
    }

    double TwoBodySE::covarianceNoCutoffs(const double &r1, const double &r2) const {
        return _energyScaleSquared * exp(-pow(r1-r2, 2.0) * _inverse2ThetaSq);
    }

    double TwoBodySE::covariance(const TwoBodyDescriptorData &r1, const TwoBodyDescriptorData &r2) {
        return covarianceNoCutoffs(r1.r, r2.r) * r1.fCut * r2.fCut;
    }

    double TwoBodySE::derivativeNoCutoffs(const double &changingR, const double &constR) const {
        return (constR - changingR) * _inverseThetaSq * covarianceNoCutoffs(changingR, constR);
    }
}
