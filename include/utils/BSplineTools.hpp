#ifndef JGAP_BSPLINETOOLS_HPP
#define JGAP_BSPLINETOOLS_HPP

#include "data/TabulationData.hpp"

namespace jgap {
    class BSplineTools {
    public:
        template<typename VectorN>
        struct InterpolationResults {
            double value;
            VectorN gradient;
        };

        static vector<double> toSplineCoefficients(const vector<double>& original, double spacing);
        static Grid1d toSplineCoefficients(const Grid1d& original);
        static Grid3d toSplineCoefficients(const Grid3d& original);

        static InterpolationResults<double> interpolate(const Grid1d& table, double pos);
        static InterpolationResults<Vector3> interpolate(const Grid3d& table, Vector3 pos);
    };
}

#endif