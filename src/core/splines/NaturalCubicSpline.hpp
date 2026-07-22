#ifndef JGAP_NATURALCUBICSPLINE_HPP
#define JGAP_NATURALCUBICSPLINE_HPP
#include <vector>

#include "Spline.hpp"
#include "core/Real.hpp"

namespace jgap {
    class NaturalCubicSpline : public Spline<1> {
    public:
        NaturalCubicSpline(const std::vector<Real>& r_vec, const std::vector<Real>& e_vec);

        InterpolationResults<1> interpolate(std::array<Real, 1> pos) const override;
        std::array<Real, 1> getCutoff() const override { return {r_vec.back()}; };

        NaturalCubicSpline* clone() const override { return new NaturalCubicSpline(*this); }

        const std::vector<Real>& getRVec() const { return r_vec; }
        const std::vector<Real>& getEnergies() const { return energies; }

    private:
        std::vector<Real> r_vec;
        std::vector<Real> energies;

        std::vector<Real> b;
        std::vector<Real> c;
        std::vector<Real> d;

        void init(const std::vector<Real>& r, const std::vector<Real>& e);
        std::size_t findInterval(Real r) const;
    };

}

#endif
