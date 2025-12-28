#ifndef JGAP_TABKERNEL3B_HPP
#define JGAP_TABKERNEL3B_HPP

#include "ThreeBodyKernel.hpp"
#include "../../../data/tabulation/Grid3d.hpp"
#include "utils/BSplineTools.hpp"
#include "utils/Utils.hpp"

namespace jgap {
    class TabKernel3b : public ThreeBodyKernel {
    public:
        TabKernel3b(const SpeciesTriplet &idSpecies, const std::shared_ptr<Grid3d> &splineGrid)
            : _idSpecies(idSpecies),
              _splineGrid(splineGrid) {
            coefficient = 1.0;
        }

        std::string getType() override { IO_NOT_INTENDED(TabKernel3b.getType); }
        DataNode serialize() override { IO_NOT_INTENDED(TabKernel3b.serialize); }
        double crossCovariance(const std::shared_ptr<IKernel> &kernel) override {
            throw std::logic_error("TabKernel3b not intended for fitting | crossCovariance()");
        }

        SpeciesTriplet getFilter() override { return _idSpecies; }

    protected:
        double valueInternal(const Vector3 &q) const override {
            return BSplineTools::interpolate(*_splineGrid, fromInvariantTriplet(q)).value;
        }
        Vector3 gradientInternal(const Vector3 &q) const override {
            Vector3 mu = fromInvariantTriplet(q);
            Vector3 dE_dmu = BSplineTools::interpolate(*_splineGrid, mu).gradient;
            std::array<Vector3, 3> _dMu_dq = dMu_dq(q);
            return _dMu_dq[0] * dE_dmu.x + _dMu_dq[1] * dE_dmu.y + _dMu_dq[2] * dE_dmu.z;
        }

    private:
        SpeciesTriplet _idSpecies;
        std::shared_ptr<Grid3d> _splineGrid;

        static Vector3 fromInvariantTriplet(const Vector3 &q) {

            const double D = std::sqrt(q.y);  // |r01 - r02|

            // two-body radii
            const double r01 = (q.x + D) * 0.5;
            const double r02 = (q.x - D) * 0.5;

            // Angle via law of cosines
            const double cosTheta = (r01*r01 + r02*r02 - q.z*q.z) / (2.0 * r01 * r02);

            return {r01, r02, cosTheta};
        }

        static std::array<Vector3, 3> dMu_dq(Vector3 q) {
            // const double r12 = q.z;

            const double D = std::sqrt(q.y);
            const double r01 = 0.5*(q.x + D);
            const double r02 = 0.5*(q.z - D);

            const double nominator = r01*r01 + r02*r02 - q.z*q.z;
            const double denominator = 2.0 * r01 * r02;

            auto dcos_dqi = [&](double dr01_dqi, double dr02_dqi, double dr12_dqi) {
                double dNominator_dqi = 2*r01*dr01_dqi + 2*r02*dr02_dqi - 2*q.z*dr12_dqi;
                double dDenominator_dqi = 2*(r02*dr01_dqi + r01*dr02_dqi);
                return (dNominator_dqi*denominator - nominator*dDenominator_dqi) / (denominator*denominator);
            };

            std::array<Vector3, 3> res{};

            // ---------------- qx ----------------
            // dr01_dqx = 0.5;
            // dr02_dqx = 0.5;
            // dr12_dqx = 0.0;
            res[0] = {0.5, 0.5, dcos_dqi(0.5, 0.5, 0.0)};

            // ---------------- qy = D² ----------------
            // dr01_dqy =  1.0 / (4.0*D);
            // dr02_dqy = -1.0 / (4.0*D);
            // dr12_dqy = 0.0
            res[1] = {1.0 / (4.0*D), -1.0 / (4.0*D),dcos_dqi(1.0 / (4.0*D), -1.0 / (4.0*D), 0.0)};

            // ---------------- qz = r12 ----------------
            // dr01_dqz = 0
            // dr02_dqz = 0
            // dr12_dqz = 1.0
            res[2] = {0.0, 0.0, dcos_dqi(0.0, 0.0, 1.0)};

            return res;
        }
    };
}

#endif
