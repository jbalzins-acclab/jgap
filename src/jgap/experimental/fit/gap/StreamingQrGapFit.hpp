#ifndef JGAP_STREAMINGQRGAPFIT_HPP
#define JGAP_STREAMINGQRGAPFIT_HPP

#include "jgap/ext/fit/gap/QRGapFit.hpp"

namespace jgap {
    class StreamingQrGapFit : public QRGapFit {
    public:
        explicit StreamingQrGapFit(Real jitter, double approx_ram_limit_gb);
        explicit StreamingQrGapFit(double approx_ram_limit_gb) : StreamingQrGapFit(1e-8, approx_ram_limit_gb) {}

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
