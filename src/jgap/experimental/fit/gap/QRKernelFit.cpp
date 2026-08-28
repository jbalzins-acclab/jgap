#include "QRKernelFit.hpp"

#include <atomic>
#include <cassert>
#include <cmath>

#include "jgap/core/UnseqFor.hpp"
#include "jgap/core/linalg/Linalg.hpp"

namespace jgap {

    std::vector<Real> QRKernelFit::findCoefficients(
        std::vector<ValuePtr<GapComponent>>& gap_components,
        const std::vector<Atoms>& training_data,
        std::vector<EnergyData>& energies_without_external,
        std::vector<Regularization>& sigmas_inverse
    ) {
        JGAP_LOG_INFO("Forming matrix A (LK_NM)");
        auto A = formMatrixA(gap_components, training_data, energies_without_external, sigmas_inverse);

        JGAP_LOG_INFO("Forming feature vector b");
        auto b = formVectorB(gap_components, energies_without_external, sigmas_inverse);

        JGAP_LOG_INFO("Solving least squares");
        auto c = leastSquares(A, b);

        return c;
    }

    std::vector<Real> QRKernelFit::leastSquares(Matrix<ColumnMajor>& A, std::vector<Real>& b) {
        return linalg::solveLeastSquaresHouseholderQR(A, b);
    }

    Matrix<ColumnMajor> QRKernelFit::formMatrixA(
        const std::vector<ValuePtr<GapComponent>>& gap_components,
        const std::vector<Atoms>& training_data,
        const std::vector<EnergyData>& energy_data,
        const std::vector<Regularization>& sigmas_inverse
    ) const {
        struct StructureData {
            const Atoms& atoms;
            const EnergyData& energy_data;
            const Regularization& sigmas_inverse;
        };

        size_t r = 0;
        std::vector<std::pair<size_t, StructureData>> starting_rows_LK_NM;
        for (size_t i = 0; i < training_data.size(); i++) {
            starting_rows_LK_NM.emplace_back(r, StructureData{training_data[i], energy_data[i], sigmas_inverse[i]});
            if (energy_data[i].energy.has_value()) r += 1;
            if (energy_data[i].forces.has_value()) r += 3 * energy_data[i].forces->size();
            if (energy_data[i].virials.has_value()) r += 6;
        }

        size_t c = 0;
        for (size_t i = 0; i < gap_components.size(); i++) {
            c += gap_components[i]->nSparsePoints();
        }

        JGAP_LOG_INFO(
            "Forming in-memory {}x{}(~{}GB) A matrix", r, c, r * c * sizeof(double) / 1024.0 / 1024.0 / 1024.0
        );
        Matrix<ColumnMajor> resulting_A(r, c);

        std::atomic counter(0);
        unseqForEach(
            starting_rows_LK_NM.begin(),
            starting_rows_LK_NM.end(),
            [&](const std::pair<size_t, StructureData>& structId) {
                auto& [starting_row, struct_data] = structId;

                size_t progress = ++counter;
                if (progress % std::max(starting_rows_LK_NM.size() / 100, 1uz) == 0) {
                    JGAP_LOG_INFO(
                        "LK_NM matrix formation progress: {} of {} ({}%)",
                        progress,
                        starting_rows_LK_NM.size(),
                        progress * 100 / starting_rows_LK_NM.size()
                    );
                }

                fillInverseSigmaLK_NM(
                    gap_components,
                    struct_data.atoms,
                    struct_data.energy_data,
                    struct_data.sigmas_inverse,
                    resulting_A,
                    starting_row
                );
            }
        );

        return resulting_A;
    }

    std::vector<Real> QRKernelFit::formVectorB(
        const std::vector<ValuePtr<GapComponent>>& components,
        const std::vector<EnergyData>& energy_data,
        const std::vector<Regularization>& sigmas_inverse
    ) {
        std::vector<Real> b;
        for (size_t i = 0; i < energy_data.size(); i++) {
            if (energy_data[i].energy.has_value()) {
                assert(sigmas_inverse[i].energy.has_value());
                b.push_back(energy_data[i].energy.value() * sigmas_inverse[i].energy.value());
            }

            if (energy_data[i].forces.has_value()) {
                assert(sigmas_inverse[i].forces.has_value());
                assert(sigmas_inverse[i].forces->size() == energy_data[i].forces->size());

                for (size_t j = 0; j < energy_data[i].forces->size(); j++) {
                    b.push_back(energy_data[i].forces->at(j).x * sigmas_inverse[i].forces->at(j).x);
                    b.push_back(energy_data[i].forces->at(j).y * sigmas_inverse[i].forces->at(j).y);
                    b.push_back(energy_data[i].forces->at(j).z * sigmas_inverse[i].forces->at(j).z);
                }
            }

            if (energy_data[i].virials.has_value()) {
                assert(sigmas_inverse[i].virials.has_value());
                b.push_back(energy_data[i].virials->xx * sigmas_inverse[i].virials->xx);
                b.push_back(energy_data[i].virials->xy * sigmas_inverse[i].virials->xy);
                b.push_back(energy_data[i].virials->xz * sigmas_inverse[i].virials->xz);
                b.push_back(energy_data[i].virials->yy * sigmas_inverse[i].virials->yy);
                b.push_back(energy_data[i].virials->yz * sigmas_inverse[i].virials->yz);
                b.push_back(energy_data[i].virials->zz * sigmas_inverse[i].virials->zz);
            }
        }

        return b;
    }

    void QRKernelFit::fillInverseSigmaLK_NM(
        const std::vector<ValuePtr<GapComponent>>& gap_components,
        const Atoms& atoms,
        const EnergyData& energy_data,
        const Regularization& sigmas_inverse,
        Matrix<ColumnMajor>& A,
        size_t starting_row
    ) {
        std::map<Real, NeighbourLists> neighbour_lists;
        for (const auto& gap_component: gap_components) {
            Real cutoff = gap_component->getCutoff();
            if (!neighbour_lists.contains(cutoff)) {
                neighbour_lists.insert({cutoff, NeighbourLists(atoms, cutoff)});
            }
        }

        size_t contribution_column = 0;
        for (const auto& gap_component: gap_components) {
            const auto& neighbour_list = neighbour_lists.at(gap_component->getCutoff());

            std::optional<AtomicQuantities> quantities = gap_component->covariate(neighbour_list);

            for (size_t point_id = 0; point_id < gap_component->nSparsePoints(); point_id++) {
                size_t row_counter = starting_row;

                if (energy_data.energy.has_value()) {
                    assert(sigmas_inverse.energy.has_value());
                    Real sigma_inverse = sigmas_inverse.energy.value();

                    Real energy = 0;
                    if (quantities.has_value()) energy = quantities->energy(point_id);

                    A(row_counter, contribution_column) = energy * sigma_inverse;
                    row_counter++;
                }

                if (energy_data.forces.has_value()) {
                    assert(sigmas_inverse.forces.has_value());
                    const auto& sigma_inverse = sigmas_inverse.forces.value();

                    for (size_t atom_id = 0; atom_id < atoms.nAtoms(); atom_id++) {
                        Vector3 force(0, 0, 0);
                        if (quantities.has_value()) force = quantities->force(point_id, atom_id);

                        A(row_counter, contribution_column) = force.x * sigma_inverse[atom_id].x;
                        row_counter++;
                        A(row_counter, contribution_column) = force.y * sigma_inverse[atom_id].y;
                        row_counter++;
                        A(row_counter, contribution_column) = force.z * sigma_inverse[atom_id].z;
                        row_counter++;
                    }
                }

                if (energy_data.virials.has_value()) {
                    assert(sigmas_inverse.virials.has_value());
                    const auto& sigma_inverse = sigmas_inverse.virials.value();

                    Virials virials(0, 0, 0, 0, 0, 0);
                    if (quantities.has_value()) virials = quantities->virials(point_id);

                    A(row_counter, contribution_column) = virials.xx * sigma_inverse.xx;
                    row_counter++;
                    A(row_counter, contribution_column) = virials.xy * sigma_inverse.xy;
                    row_counter++;
                    A(row_counter, contribution_column) = virials.xz * sigma_inverse.xz;
                    row_counter++;
                    A(row_counter, contribution_column) = virials.yy * sigma_inverse.yy;
                    row_counter++;
                    A(row_counter, contribution_column) = virials.yz * sigma_inverse.yz;
                    row_counter++;
                    A(row_counter, contribution_column) = virials.zz * sigma_inverse.zz;
                    row_counter++;
                }

                contribution_column++;
            }
        }
    }
}
