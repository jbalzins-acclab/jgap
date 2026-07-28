#include "CubicBSpline3D.hpp"
#include "CubicBSpline.hpp"

#include <cmath>

namespace jgap {

    CubicBSpline3D::CubicBSpline3D(const Grid<3>& coefficients) : coefficients(coefficients) {}

    CubicBSpline3D CubicBSpline3D::fit(const Grid<3>& values) {
        const size_t Mx = values.sizes[0];
        const size_t My = values.sizes[1];
        const size_t Mz = values.sizes[2];

        const size_t Nx = Mx + 2;
        const size_t Ny = My + 2;
        const size_t Nz = Mz + 2;

        std::array<Real, 3> new_origin = values.origin;
        for (size_t i = 0; i < 3; i++) {
            new_origin[i] -= values.spacing[i];
        }

        Grid<3> temp1({Nx, My, Mz}, values.spacing, new_origin);
        for (size_t iy = 0; iy < My; ++iy) {
            for (size_t iz = 0; iz < Mz; ++iz) {
                std::vector<Real> slice(Mx);
                for (size_t ix = 0; ix < Mx; ++ix) {
                    slice[ix] = values({ix, iy, iz});
                }
                std::vector<Real> coeffs = CubicBSpline::toSplineCoefficients(slice, values.spacing[0]);
                for (size_t ix = 0; ix < Nx; ++ix) {
                    temp1({ix, iy, iz}) = coeffs[ix];
                }
            }
        }

        Grid<3> temp2({Nx, Ny, Mz}, values.spacing, new_origin);
        for (size_t ix = 0; ix < Nx; ++ix) {
            for (size_t iz = 0; iz < Mz; ++iz) {
                std::vector<Real> slice(My);
                for (size_t iy = 0; iy < My; ++iy) {
                    slice[iy] = temp1({ix, iy, iz});
                }
                std::vector<Real> coeffs = CubicBSpline::toSplineCoefficients(slice, values.spacing[1]);
                for (size_t iy = 0; iy < Ny; ++iy) {
                    temp2({ix, iy, iz}) = coeffs[iy];
                }
            }
        }

        Grid<3> final_coeff({Nx, Ny, Nz}, values.spacing, new_origin);
        for (size_t ix = 0; ix < Nx; ++ix) {
            for (size_t iy = 0; iy < Ny; ++iy) {
                std::vector<Real> slice(Mz);
                for (size_t iz = 0; iz < Mz; ++iz) {
                    slice[iz] = temp2({ix, iy, iz});
                }
                std::vector<Real> coeffs = CubicBSpline::toSplineCoefficients(slice, values.spacing[2]);
                for (size_t iz = 0; iz < Nz; ++iz) {
                    final_coeff({ix, iy, iz}) = coeffs[iz];
                }
            }
        }

        return CubicBSpline3D(final_coeff);
    }

    InterpolationResults<3> CubicBSpline3D::interpolate(std::array<Real, 3> pos) const {
        const auto& coeffs = coefficients.data_flat;
        const auto& dims = coefficients.sizes;
        const auto& spacing = coefficients.spacing;
        const auto& origin = coefficients.origin;
        const auto cutoff = getCutoff();

        for (size_t d = 0; d < 3; d++) {
            const Real data_origin = origin[d] + spacing[d];
            const Real data_upper = data_origin + static_cast<Real>(dims[d] - 3) * spacing[d];
            if (pos[d] < (data_origin - 1e-9) || pos[d] > (data_upper + 1e-9)) {
                return {0.0, {0.0, 0.0, 0.0}};
            }
        }

        Real u[3], t[3];
        int ii[3];

        for (size_t d = 0; d < 3; d++) {
            const Real data_origin = origin[d] + spacing[d];
            u[d] = (pos[d] - data_origin) / spacing[d];
            int i0 = static_cast<int>(std::floor(u[d]));
            int imax = static_cast<int>(dims[d]) - 4;
            ii[d] = std::clamp(i0, 0, imax);
            t[d] = u[d] - ii[d];
        }

        Real Phi[3][4], dPhi[3][4];
        for (size_t d = 0; d < 3; ++d) {
            const Real t2 = t[d] * t[d];
            const Real t3 = t2 * t[d];
            const Real dinv = 1.0 / spacing[d];

            Phi[d][0] = (1 - 3 * t[d] + 3 * t2 - t3) / 6.0;
            Phi[d][1] = (4 - 6 * t2 + 3 * t3) / 6.0;
            Phi[d][2] = (1 + 3 * t[d] + 3 * t2 - 3 * t3) / 6.0;
            Phi[d][3] = t3 / 6.0;

            dPhi[d][0] = (-3 + 6 * t[d] - 3 * t2) * dinv / 6.0;
            dPhi[d][1] = (-12 * t[d] + 9 * t2) * dinv / 6.0;
            dPhi[d][2] = (3 + 6 * t[d] - 9 * t2) * dinv / 6.0;
            dPhi[d][3] = (3 * t2) * dinv / 6.0;
        }

        const int N1 = dims[1];
        const int N2 = dims[2];

        Real value = 0.0;
        Real dval[3] = {0.0, 0.0, 0.0};

        for (int i = 0; i < 4; ++i) {
            const int base_i = ((ii[0] + i) * N1 + ii[1]) * N2 + ii[2];

            Real ppc = 0.0;
            Real dpp1 = 0.0;
            Real dpp2 = 0.0;

            for (int j = 0; j < 4; ++j) {
                const Real* cptr = &coeffs[base_i + j * N2];

                const Real pc = Phi[2][0] * cptr[0] + Phi[2][1] * cptr[1] + Phi[2][2] * cptr[2] + Phi[2][3] * cptr[3];
                const Real dpc =
                    dPhi[2][0] * cptr[0] + dPhi[2][1] * cptr[1] + dPhi[2][2] * cptr[2] + dPhi[2][3] * cptr[3];

                ppc += Phi[1][j] * pc;
                dpp1 += dPhi[1][j] * pc;
                dpp2 += Phi[1][j] * dpc;
            }

            value += Phi[0][i] * ppc;
            dval[0] += dPhi[0][i] * ppc;
            dval[1] += Phi[0][i] * dpp1;
            dval[2] += Phi[0][i] * dpp2;
        }

        return {value, {dval[0], dval[1], dval[2]}};
    }

    std::array<Real, 3> CubicBSpline3D::getCutoff() const {
        // Original grid data_origin = origin + spacing, data_points N = sizes - 2.
        // Cutoff is exact upper bound: data_origin + (N - 1) * spacing.
        std::array<Real, 3> cutoff;
        for (size_t d = 0; d < 3; d++) {
            const Real data_points = static_cast<Real>(coefficients.sizes[d] - 2);
            const Real data_origin = coefficients.origin[d] + coefficients.spacing[d];
            cutoff[d] = data_origin + (data_points - 1) * coefficients.spacing[d];
        }
        return cutoff;
    }
}
