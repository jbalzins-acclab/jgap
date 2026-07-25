#ifndef JGAP_STREAMINGQRGAPFIT_HPP
#define JGAP_STREAMINGQRGAPFIT_HPP

#include "jgap/ext/fit/gap/QRGapFit.hpp"

namespace jgap {
    class StreamingQrGapFit : public QRGapFit {
    public:
        explicit StreamingQrGapFit(
            Real jitter = 1e-8,
            size_t target_chunk_rows = 8000
        );

    protected:
        std::vector<Real> findCoefficients(
            std::vector<ValuePtr<GapComponent>>& gap_components,
            const std::vector<Atoms>& training_data,
            std::vector<EnergyData>& energies_without_external,
            std::vector<Regularization>& sigmas_inverse
        ) override;

    protected:
        size_t target_chunk_rows{8000};
    };
}

#endif
