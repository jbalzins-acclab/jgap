#ifndef JGAP_BLOCKINCREMENTALQRGAPFIT_HPP
#define JGAP_BLOCKINCREMENTALQRGAPFIT_HPP

#include "jgap/ext/fit/gap/QRGapFit.hpp"

namespace jgap {
    class BlockIncrementalQRGapFit : public QRGapFit {
    public:
        explicit BlockIncrementalQRGapFit(Real jitter, double approx_ram_limit_gb);
        explicit BlockIncrementalQRGapFit(double approx_ram_limit_gb) : BlockIncrementalQRGapFit(1e-8, approx_ram_limit_gb) {}

    protected:
        std::vector<Real> findCoefficients(
            std::vector<ValuePtr<GapComponent>>& gap_components,
            const std::vector<Atoms>& training_data,
            std::vector<EnergyData>& energies_without_external,
            std::vector<Regularization>& sigmas_inverse
        ) override;

    protected:
        double approx_ram_limit_gb;
    };
}

#endif
