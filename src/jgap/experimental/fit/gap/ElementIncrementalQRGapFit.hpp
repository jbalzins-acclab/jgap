#ifndef JGAP_ELEMENTINCREMENTALQRGAPFIT_HPP
#define JGAP_ELEMENTINCREMENTALQRGAPFIT_HPP

#include <memory>
#include <optional>
#include <set>
#include <vector>

#include "jgap/core/atomic/species/Species.hpp"
#include "jgap/core/linalg/BlockIncrementalQRAccumulator.hpp"
#include "jgap/ext/fit/gap/QRGapFit.hpp"

namespace jgap {

    class ElementIncrementalQRGapFit : public QRGapFit {
    public:
        explicit ElementIncrementalQRGapFit(double jitter, double approx_ram_limit_gb);

        explicit ElementIncrementalQRGapFit(double approx_ram_limit_gb) :
            ElementIncrementalQRGapFit(1e-8, approx_ram_limit_gb) {}

    protected:
        void findCoefficients(
            GapPotential& to_be_fit,
            const std::vector<Atoms>& training_data,
            std::vector<EnergyData>& energies_without_external,
            std::vector<Regularization>& sigmas_inverse
        ) override;

    private:
        struct IncrementMeta {
            Species species = Species::Anon();
            std::vector<ValuePtr<GapComponent>> single_element_components;
            std::vector<size_t> single_element_structure_indices;

            std::vector<ValuePtr<GapComponent>> new_multi_element_components;
            std::vector<size_t> new_multi_element_structure_indices;
        };

        double approx_ram_limit_gb;

        SharedArray<double> allocateMatrixMemory(const std::vector<ValuePtr<GapComponent>>& gap_components) const;

        std::vector<IncrementMeta> determineIncrements(
            std::vector<ValuePtr<GapComponent>>& gap_components,
            const std::vector<Atoms>& training_data,
            std::vector<size_t>& original_component_indices
        );

        std::vector<double> doQRMakeContagiousAndReturnAccumulatedB(
            Matrix<ColumnMajor> A,
            SharedArray<double> b,
            size_t n_cols,
            const std::vector<ValuePtr<GapComponent>>& components,
            const std::vector<size_t>& structure_indices,
            const std::vector<Atoms>& training_data,
            const std::vector<EnergyData>& energy_data,
            const std::vector<Regularization>& sigmas_inverse
        );

        std::vector<double> singleElementQR(
            SharedArray<double>& single_element_memory_space,
            const std::vector<ValuePtr<GapComponent>>& components,
            const std::vector<size_t>& single_element_structure_indices,
            const std::vector<Atoms>& training_data,
            const std::vector<EnergyData>& energy_data,
            const std::vector<Regularization>& sigmas_inverse
        );

        Matrix<ColumnMajor> restackContagiousBlocks(
            SharedArray<double>& full_memory, size_t former_M, size_t single_element_cols, size_t require_n_cols
        );

        void addUMMs(
            Matrix<ColumnMajor>& A,
            size_t filled_M,
            const std::vector<ValuePtr<GapComponent>>& new_multi_element_components
        );
    };
}

#endif
