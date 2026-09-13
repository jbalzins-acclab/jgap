#ifndef JGAP_NATURALCUBICSPLINE_HPP
#define JGAP_NATURALCUBICSPLINE_HPP
#include <vector>

#include "Spline.hpp"

namespace jgap {
    class NaturalCubicSpline : public Spline<1> {
    public:
        NaturalCubicSpline(const std::vector<double>& r_vec, const std::vector<double>& e_vec);

        InterpolationResults<1> interpolate(std::array<double, 1> pos) const override;
        std::array<double, 1> getCutoff() const override { return {r_vec.back()}; };

        NaturalCubicSpline* clone() const override { return new NaturalCubicSpline(*this); }

        const std::vector<double>& getRVec() const { return r_vec; }
        const std::vector<double>& getEnergies() const { return energies; }

    private:
        std::vector<double> r_vec;
        std::vector<double> energies;

        std::vector<double> b;
        std::vector<double> c;
        std::vector<double> d;

        void init(const std::vector<double>& r, const std::vector<double>& e);
        std::size_t findInterval(double r) const;
    };

}

#endif
