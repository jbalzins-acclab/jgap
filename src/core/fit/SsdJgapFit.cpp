#ifndef _WIN32
#include "core/fit/SsdJgapFit.hpp"

#include "core/matrices/sigmas/SimpleRegularizationRules.hpp"
#include "core/neighbours/NeighbourFinder.hpp"
#include "io/log/StdoutLogger.hpp"

#include <tbb/parallel_for_each.h>

#include "utils/Utils.hpp"

namespace jgap {

    SsdJgapFit::SsdJgapFit(const nlohmann::json &params) {

        _descriptors = {};
        for (const auto& [label, descriptor] : params["descriptors"].items()) {
            _descriptors[label] = (ParserRegistry<Descriptor>::get(descriptor));
        }

        _regularizationRules = ParserRegistry<RegularizationRules>::get(params["sigma_rules"]);
        _jitter = params.value("jitter", 1e-8);
    }

    shared_ptr<Potential> SsdJgapFit::fit(const vector<AtomicStructure>& trainingData) {
        CurrentLogger::get()->info("Starting JGAP fit");

        vector _trainingData(trainingData);

        double maxCutoff = 0;

        CurrentLogger::get()->info("Checking sparse points");

        auto descriptorsAsVec = vector<shared_ptr<Descriptor>>();
        for (const auto& descriptor: _descriptors | views::values) {
            descriptorsAsVec.push_back(descriptor);

            // To avoid ugly "cutoff" in sparse json
            NeighbourFinder::findNeighbours(_trainingData, descriptor->getCutoff());
            if (descriptor->nSparsePoints() == 0) {
                CurrentLogger::get()->info(format(
                    "No sparse points found in descriptor {} -> setting from data",
                    descriptor->serialize().dump()
                    ));
                descriptor->setSparsePoints(_trainingData);
            }

            maxCutoff = max(maxCutoff, descriptor->getCutoff());
        }
        CurrentLogger::get()->info("Sparsification complete");

        CurrentLogger::get()->debug("Full neighbour-list");
        NeighbourFinder::findNeighbours(_trainingData, maxCutoff);

        //// ----------------------------------------------------------------------------------------------------

        CurrentLogger::get()->info("Per-structure regularization setup");
        for (auto& structure: _trainingData) {
            _regularizationRules->fillSigmas(structure);
        }
        CurrentLogger::get()->info("Per-structure regularization setup finished");

        //// ----------------------------------------------------------------------------------------------------

        CurrentLogger::get()->info("Making matrix A");
        auto A = makeA(descriptorsAsVec, _trainingData);
        CurrentLogger::get()->info("Done making matrix A");
        // CurrentLogger::get()->info(matrixToString(A));

        CurrentLogger::get()->info("Making feature vector b");
        auto b = makeB(descriptorsAsVec, _trainingData);
        CurrentLogger::get()->info("Done making feature vector b");

        //// ----------------------------------------------------------------------------------------------------

        CurrentLogger::get()->info("Doing linear algebra");
        vector c = leastSquares(A, b);
        CurrentLogger::get()->info("Finished linear algebra");

        size_t counter = 0;
        for (const auto &descriptor: _descriptors | views::values) {
            const size_t n = descriptor->nSparsePoints();

            // auto bb = vector<double>{c.data() + counter, c.data() + counter + n};
            descriptor->setCoefficients(vector<double>(c.begin() + counter, c.begin() + counter + n));

            counter += n;
        }

        return make_shared<JgapPotential>(_descriptors);
    }

    vector<double> SsdJgapFit::leastSquares(Eigen::MatrixXd &A, Eigen::VectorXd &b) {
        CurrentLogger::get()->debug("Init Eigen::HouseholderQR");
        const Eigen::HouseholderQR<Eigen::Ref<Eigen::MatrixXd>> qr(A);

        CurrentLogger::get()->debug("Q^t");
        auto Qt = qr.householderQ().transpose();

        CurrentLogger::get()->debug("Q^t * b");
        auto Qt_b =  Qt * b;

        CurrentLogger::get()->debug("R");
        auto R = qr.matrixQR().topLeftCorner(A.cols(), A.cols());

        CurrentLogger::get()->debug("R^-1 * Qt_b");
        Eigen::VectorXd c = R.triangularView<Eigen::Upper>().solve(Qt_b.head(A.cols()));
        CurrentLogger::get()->debug("c" + to_string(c[0]));

        return vector<double>{c.data(), c.data() + c.size()};
    }

    Eigen::MatrixXd SsdJgapFit::makeA(const vector<shared_ptr<Descriptor>> &descriptors,
                                        const vector<AtomicStructure> &atomicStructures) const {

        size_t r = 0;
        vector<pair<size_t, AtomicStructure>> startingRowsK_nm;
        for (const auto& structure : atomicStructures) {
            startingRowsK_nm.emplace_back(r, structure);
            if (structure.energy.has_value()) r += 1;
            if (structure.forces.has_value()) r += 3 * structure.size();
            if (structure.virials.has_value()) r += 6;
        }

        vector<tuple<size_t/*row*/, size_t/*col*/, size_t/*desc_idx*/>> startingPointsK_mm;
        size_t c = 0;
        for (size_t i = 0; i < descriptors.size(); i++) {
            startingPointsK_mm.emplace_back(r + c, c, i);
            c += descriptors[i]->nSparsePoints();
        }

        CurrentLogger::get()->info("Forming in-memory " + to_string(r+c) + "x" + to_string(c) + " A matrix");
        Eigen::MatrixXd resultingA = Eigen::MatrixXd::Zero(r + c, c);

        atomic counter(0);
        tbb::parallel_for_each(
            startingRowsK_nm.begin(), startingRowsK_nm.end(),
            [&](const pair<size_t, AtomicStructure>& structId) {

                size_t progress = ++counter;
                if (progress % max(startingRowsK_nm.size() / 100, 1uz) == 0) {
                    CurrentLogger::get()->debug(format(
                            "K_nm matrix formation progress: {} of {} ({}%)",
                            progress,
                            startingRowsK_nm.size(),
                            progress * 100 / startingRowsK_nm.size()
                            ));
                }

                fillInverseSigmaK_nm(descriptors, structId.second, resultingA, structId.first);
            }
        );

        tbb::parallel_for_each(
            startingPointsK_mm.begin(), startingPointsK_mm.end(),
            [&](const tuple<size_t/*row*/, size_t/*col*/, size_t/*desc_idx*/>& descriptorId) {
                CurrentLogger::get()->debug(format("K_mm for descriptor {}", get<2>(descriptorId)));
                fillU_mm(get<0>(descriptorId), get<1>(descriptorId), *descriptors[get<2>(descriptorId)], resultingA);
            }
        );

        return resultingA;
    }

    Eigen::VectorXd SsdJgapFit::makeB(const vector<shared_ptr<Descriptor>> &descriptors,
                                                    const vector<AtomicStructure> &atomicStructures) {
        vector<double> b;
        for (auto& structure: atomicStructures) {
            if (structure.energy.has_value()) {
                b.push_back(structure.energy.value() * structure.energySigmaInverse.value());
            }
            if (structure.forces.has_value()) {
                for (const auto& atom: structure) {
                    b.push_back(atom.force().x * atom.forceSigmasInverse().x);
                    b.push_back(atom.force().y * atom.forceSigmasInverse().y);
                    b.push_back(atom.force().z * atom.forceSigmasInverse().z);
                }
            }
            if (structure.virials.has_value()) {
                b.push_back(structure.virials.value()[0].x * structure.virialSigmasInverse.value()[0].x);
                b.push_back(structure.virials.value()[0].y * structure.virialSigmasInverse.value()[0].y);
                b.push_back(structure.virials.value()[0].z * structure.virialSigmasInverse.value()[0].z);

                b.push_back(structure.virials.value()[1].y * structure.virialSigmasInverse.value()[1].y);
                b.push_back(structure.virials.value()[1].z * structure.virialSigmasInverse.value()[1].z);

                b.push_back(structure.virials.value()[2].z * structure.virialSigmasInverse.value()[2].z);
            }
        }

        for (auto& descriptor: descriptors) {
            b.resize(b.size() + descriptor->nSparsePoints());
        }

        return Eigen::Map<Eigen::VectorXd>(b.data(), b.size());
    }

    void SsdJgapFit::fillInverseSigmaK_nm(const vector<shared_ptr<Descriptor>> &descriptors,
                                                const AtomicStructure &atomicStructure,
                                                Eigen::MatrixXd &A,
                                                const size_t startingRow) {
        size_t contributionColumn = 0;
        for (const auto& descriptor : descriptors) {
            auto contributions = descriptor->covariate(atomicStructure);

            for (const auto& contribution: contributions) {

                size_t currentRow = startingRow;

                if (atomicStructure.energy.has_value()) {
                    A(currentRow++, contributionColumn) = contribution.total
                                                                    * atomicStructure.energySigmaInverse.value();
                }

                if (atomicStructure.forces.has_value()) {
                    for (size_t rowInc = 0; rowInc < contribution.forces.size(); rowInc++) {

                        const auto force = contribution.forces[rowInc]; // FORCE IS NEGATIVE
                        const auto fSigmasInverse = (*atomicStructure.forceSigmasInverse)[rowInc];

                        A(currentRow++, contributionColumn) = force.x * fSigmasInverse.x;
                        A(currentRow++, contributionColumn) = force.y * fSigmasInverse.y;
                        A(currentRow++, contributionColumn) = force.z * fSigmasInverse.z;
                    }
                }

                if (atomicStructure.virials.has_value()) {
                    array virials = contribution.virials;
                    A(currentRow++, contributionColumn) = virials[0].x
                        * atomicStructure.virialSigmasInverse.value()[0].x;
                    A(currentRow++, contributionColumn) = virials[0].y
                        * atomicStructure.virialSigmasInverse.value()[0].y;
                    A(currentRow++, contributionColumn) = virials[0].z
                        * atomicStructure.virialSigmasInverse.value()[0].z;

                    A(currentRow++, contributionColumn) = virials[1].y
                        * atomicStructure.virialSigmasInverse.value()[1].y;
                    A(currentRow++, contributionColumn) = virials[1].z
                        * atomicStructure.virialSigmasInverse.value()[1].z;

                    A(currentRow++, contributionColumn) = virials[2].z
                        * atomicStructure.virialSigmasInverse.value()[2].z;
                }

                contributionColumn++;
            }
        }
    }

    void SsdJgapFit::fillU_mm(const size_t startingRow, const size_t startingCol,
                                Descriptor &descriptor, Eigen::MatrixXd &A) const {

        for (auto &[sparseId, K_mmPart]: descriptor.selfCovariate()) {

            size_t n = K_mmPart->rows();
            for (size_t i = 0; i < n; i++) (*K_mmPart)(i, i) += _jitter;

            auto K_mmConverted = convertToEigen(*K_mmPart);
            auto U_mmPart = choleskyDecomposition(K_mmConverted);

            for (size_t i = 0; i < n; i++) {
                for (size_t j = 0; j < n; j++) {
                    A(startingRow + sparseId + i, startingCol + sparseId + j) = U_mmPart(i, j);
                }
            }
        }
    }

    Eigen::MatrixXd SsdJgapFit::choleskyDecomposition(Eigen::MatrixXd &matrix) {
        Eigen::LLT<Eigen::MatrixXd> llt(matrix);
        if (llt.info() != Eigen::Success)
            CurrentLogger::get()->error("Cholesky decomposition failed: matrix not positive definite", true);

        return llt.matrixU();
    }

    Eigen::MatrixXd SsdJgapFit::convertToEigen(MatrixBlock &matrixBlock) {
        return Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>(
            matrixBlock.rawData().data(), matrixBlock.rows(), matrixBlock.columns()
            );
    }
}

#include <iostream>
#include <vector>
#include <memory>
#include <cstring>
#include <algorithm>
#include <cmath>

#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>

#include <Eigen/Dense>
#include <Eigen/QR>

using Matrix = Eigen::MatrixXd;
using Vector = Eigen::VectorXd;

class MemoryMappedMatrix {
private:
    std::string filename_;
    size_t rows_, cols_;
    size_t file_size_;

    int fd_;
    void* mapped_ptr_;

public:
    MemoryMappedMatrix(const std::string& filename, size_t rows, size_t cols, bool create = true)
        : filename_(filename), rows_(rows), cols_(cols),
          file_size_(rows * cols * sizeof(double)) {

        if (create) {
            create_and_map();
        } else {
            open_and_map();
        }
    }

    ~MemoryMappedMatrix() {
        unmap();
    }

    double* data() { return static_cast<double*>(mapped_ptr_); }

    double& operator()(size_t i, size_t j) {
        return static_cast<double*>(mapped_ptr_)[i * cols_ + j];
    }

    const double& operator()(size_t i, size_t j) const {
        return static_cast<double*>(mapped_ptr_)[i * cols_ + j];
    }

    double* block(size_t row_start, size_t col_start) {
        return &((*this)(row_start, col_start));
    }

    size_t rows() const { return rows_; }
    size_t cols() const { return cols_; }
    size_t stride() const { return cols_; }

    void sync() {
        msync(mapped_ptr_, file_size_, MS_SYNC);
    }


private:
    void create_and_map() {
        fd_ = open(filename_.c_str(), O_CREAT | O_RDWR | O_TRUNC, 0644);
        if (fd_ == -1) throw std::runtime_error("Failed to create file");

        if (ftruncate(fd_, file_size_) == -1) {
            close(fd_); throw std::runtime_error("Failed to set file size");
        }

        mapped_ptr_ = mmap(nullptr, file_size_, PROT_READ | PROT_WRITE, MAP_SHARED, fd_, 0);
        if (mapped_ptr_ == MAP_FAILED) {
            close(fd_); throw std::runtime_error("Failed to mmap file");
        }

        madvise(mapped_ptr_, file_size_, MADV_SEQUENTIAL);
        std::memset(mapped_ptr_, 0, file_size_);
    }

    void open_and_map() {
        fd_ = open(filename_.c_str(), O_RDWR);
        if (fd_ == -1) throw std::runtime_error("Failed to open file");

        mapped_ptr_ = mmap(nullptr, file_size_, PROT_READ | PROT_WRITE, MAP_SHARED, fd_, 0);
        if (mapped_ptr_ == MAP_FAILED) {
            close(fd_); throw std::runtime_error("Failed to mmap file");
        }

        madvise(mapped_ptr_, file_size_, MADV_RANDOM);
    }

    void unmap() {
        if (mapped_ptr_ != MAP_FAILED) { munmap(mapped_ptr_, file_size_); mapped_ptr_ = nullptr; }
        if (fd_ != -1) { close(fd_); fd_ = -1; }
    }
};

class OutOfCoreLeastSquares {
private:
    size_t total_rows_;  // n + m
    size_t n_;           // number of data points
    size_t m_;           // number of parameters (columns)
    size_t available_ram_doubles_;
    size_t mb_, nb_;     // Block sizes

    // std::unique_ptr<MemoryMappedMatrix> matrix_;
    std::string temp_dir_;
    std::string matrix_filename_;

    // Store Q transformations compactly
    struct QTransformation {
        size_t k;                    // Column index
        size_t row_start, row_end;   // Row range this transformation affects
        Matrix Q_panel;              // Panel transformation matrix
        Vector qtb_contribution;     // Q^T * b contribution for this block
    };
    std::vector<QTransformation> q_transformations_;

    // Final R matrix (upper triangular, fits in memory)
    Matrix R_;

public:
    OutOfCoreLeastSquares(size_t n, size_t m, size_t available_ram_doubles,
                         const std::string& temp_dir = "/tmp")
        : total_rows_(n + m), n_(n), m_(m),
          available_ram_doubles_(available_ram_doubles), temp_dir_(temp_dir) {

        calculate_block_sizes();

        matrix_filename_ = temp_dir_ + "/ols_matrix_" +
                          std::to_string(reinterpret_cast<uintptr_t>(this)) + ".dat";
        // matrix_ = std::make_unique<MemoryMappedMatrix>(matrix_filename_, total_rows_, m_, true);

        R_ = Matrix::Zero(m_, m_);  // R is m x m

        std::cout << "Initialized OOC Least Squares: (" << n_ << "+" << m_
                  << ")x" << m_ << ", blocks: " << mb_ << "x" << nb_ << std::endl;
    }

    ~OutOfCoreLeastSquares() {
        std::remove(matrix_filename_.c_str());
    }

    // Main method: solve least squares with structured RHS
    Vector solve_least_squares(shared_ptr<MemoryMappedMatrix>& A, const Vector& b_data) {
        if (b_data.size() != n_) {
            throw std::invalid_argument("b_data must have n elements");
        }

        // Create full b vector: [b_data; zeros]
        Vector b_full = Vector::Zero(total_rows_);
        b_full.head(n_) = b_data;

        // Perform blocked QR and accumulate Q^T * b
        Vector qtb = compute_qr_and_qtb(A, b_full);

        // Solve R * x = Q^T * b using back substitution
        return solve_triangular(qtb);
    }

private:
    void calculate_block_sizes() {
        size_t usable_ram = static_cast<size_t>(available_ram_doubles_ * 0.8);

        nb_ = std::min(static_cast<size_t>(64), m_);

        // Need space for: panel (mb x nb) + qtb computation + workspace
        // Approximately: mb * (nb + m) < usable_ram
        if (nb_ + m_ > 0) {
            mb_ = usable_ram / (nb_ + m_ + nb_);  // Extra nb for workspace
            mb_ = std::max(mb_, static_cast<size_t>(1));
            mb_ = std::min(mb_, total_rows_);
        } else {
            mb_ = std::min(usable_ram / 1000, total_rows_);
        }

        if (mb_ < 32 && nb_ > 8) {
            nb_ = std::max(static_cast<size_t>(8), nb_ / 2);
            mb_ = usable_ram / (nb_ + m_ + nb_);
            mb_ = std::max(mb_, static_cast<size_t>(1));
        }
    }

    Vector compute_qr_and_qtb(shared_ptr<MemoryMappedMatrix>& A, const Vector& b_full) {
        q_transformations_.clear();
        Vector qtb = b_full;  // Will be transformed in place

        for (size_t k = 0; k < m_; k += nb_) {
            size_t panel_width = std::min(nb_, m_ - k);
            if (panel_width == 0) break;

            std::cout << "Processing panel at column " << k
                      << ", width " << panel_width << std::endl;

            // Process panel in row blocks
            for (size_t row_start = k; row_start < total_rows_; row_start += mb_) {
                size_t row_end = std::min(row_start + mb_, total_rows_);
                size_t block_height = row_end - row_start;

                if (block_height == 0) break;

                // Extract panel block
                Eigen::Map<Matrix> panel_block(
                    A->block(row_start, k),
                    block_height, panel_width,
                    A->stride()
                );

                if (row_start == k) {
                    // First block - perform QR
                    Eigen::HouseholderQR<Matrix> qr(panel_block);
                    Matrix Q_panel = qr.householderQ();
                    Matrix R_panel = qr.matrixQR().triangularView<Eigen::Upper>();

                    // Store R factor
                    size_t r_rows = std::min(panel_width, static_cast<size_t>(R_panel.rows()));
                    R_.block(k, k, r_rows, panel_width) = R_panel.topRows(r_rows);

                    // Apply Q^T to corresponding part of b
                    Vector b_block = qtb.segment(row_start, block_height);
                    Vector qtb_block = Q_panel.transpose() * b_block;
                    qtb.segment(row_start, block_height) = qtb_block;

                    // Store transformation for trailing matrix updates
                    QTransformation qtrans;
                    qtrans.k = k;
                    qtrans.row_start = row_start;
                    qtrans.row_end = row_end;
                    qtrans.Q_panel = Q_panel;
                    qtrans.qtb_contribution = qtb_block;

                    q_transformations_.push_back(std::move(qtrans));

                    // Apply to trailing matrix
                    if (k + panel_width < m_) {
                        apply_panel_to_trailing(Q_panel, k, panel_width, row_start, row_end);
                    }
                } else {
                    // Subsequent blocks - apply existing Q transformations
                    // This would involve applying previous Q transformations
                    // For simplicity, we'll handle this case similarly
                    // In practice, you'd accumulate the transformations
                }
            }
        }

        return qtb;
    }

    void apply_panel_to_trailing(const Matrix& Q_panel, size_t k, size_t panel_width,
                                size_t row_start, size_t row_end) {
        size_t remaining_cols = m_ - k - panel_width;

        for (size_t col_offset = 0; col_offset < remaining_cols; col_offset += nb_) {
            size_t trailing_cols = std::min(nb_, remaining_cols - col_offset);
            size_t trailing_start_col = k + panel_width + col_offset;

            Eigen::Map<Matrix> trailing_block(
                matrix_->block(row_start, trailing_start_col),
                row_end - row_start, trailing_cols,
                matrix_->stride()
            );

            trailing_block = Q_panel.transpose() * trailing_block;
        }
    }

    Vector solve_triangular(const Vector& qtb) {
        // Solve R * x = (Q^T * b)[1:m] using back substitution
        // Only need first m elements of Q^T * b
        Vector rhs = qtb.head(m_);
        Vector x = Vector::Zero(m_);

        // Back substitution
        for (int i = m_ - 1; i >= 0; --i) {
            double sum = rhs(i);
            for (size_t j = i + 1; j < m_; ++j) {
                sum -= R_(i, j) * x(j);
            }
            x(i) = sum / R_(i, i);
        }

        return x;
    }

public:

    // Get the R matrix (for analysis)
    const Matrix& get_R() const { return R_; }
};

// Test function
void test_least_squares() {
    const size_t n = 1000;  // Data points
    const size_t m = 50;    // Parameters
    const size_t available_ram = 100000;  // 100K doubles

    Vector b_data;

    // Solve using out-of-core method
    OutOfCoreLeastSquares solver(n, m, available_ram);

    Vector x_computed = solver.solve_least_squares(A, b_data);
}

#endif
