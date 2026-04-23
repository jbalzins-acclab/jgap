#ifndef JGAP_BSPLINETOOLS_HPP
#define JGAP_BSPLINETOOLS_HPP

#include "../data/tabulation/TabulationData.hpp"

namespace jgap {
    class BSplineTools {
    public:
        template<typename VectorN>
        struct InterpolationResults {
            double value;
            VectorN gradient;
        };

        static std::vector<double> toSplineCoefficients(const std::vector<double>& original, double spacing);
        static Grid1d toSplineCoefficients(const Grid1d& original);
        static Grid3d toSplineCoefficients(const Grid3d& original);

        static InterpolationResults<double> interpolate(const Grid1d& table, double pos);
        static InterpolationResults<Vector3> interpolate(const Grid3d& table, const Vector3 &pos);
    };
}

#endif