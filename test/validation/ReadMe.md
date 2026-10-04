# JGAP Automated Validation Suite

Comprehensive automated validation benchmark suite for `jgap`.

This suite rigorously validates every stage of potential development against industry-standard reference baselines (**Fortran QUIP `gap_fit`**, **Python `tabGAP`**, **LAMMPS**, and **DFT literature values**) without requiring external binaries to be installed at test time.

---

## 1. Pipeline Overview

The validation suite validates each phase of the machine learning potential pipeline:

1. **Reference Data & Training Datasets**:
   - `test/resources/validation.zip`: Contains reference QUIP XML potentials, test set predictions (`pred.xyz`), tabGAP spline tables, LAMMPS predictions, and DFT benchmarks.
   - `test/structure-databases/`: Training databases for elemental and alloy potentials (`Al`, `Cu`, `Ni`, `Fe`, `FeNi`).

2. **Validation Benchmarks**:
   - **`ValidateQRGapFit`**: QR regression fit against reference QUIP XML coefficients.
   - **`ValidateGapEnergyEval`**: In-engine energy and force predictions against reference QUIP predictions.
   - **`ValidateTabulation`**: JGAP spline tabulation against reference tabGAP tables.
   - **`ValidateTabGapEnergyEval`**: JGAP `TabGapPotential` engine predictions against LAMMPS reference predictions.
   - **`ValidateQrVariants`**: Equivalence and performance of QR solver variants (`FullQR`, `BlockIncrementalQR`, `ElementIncrementalQR`, out-of-core streaming).
   - **`test_elemental_properties.py`**: Potential refitting via `standard_gap_fit`, property evaluations via ASE + `elastic`, and comparison against DFT literature values.

3. **Reports & Metrics**:
   - Automated execution and regression enforcement via `ctest -L validation`.
   - Metric outputs analyzed and formatted via `analyze_validation.py`.

---

## 2. Prerequisites & Setup

### 2.1 Extract Reference Validation Data
All reference files (QUIP XMLs, predictions, tabGAP tables, LAMMPS references) are archived in `test/resources/validation.zip` (52 MB).

Extract before running validation benchmarks:
```bash
unzip -q test/resources/validation.zip -d test/resources/
```

### 2.2 CMake Configuration
Ensure CMake is configured with `JGAP_BUILD_VALIDATION=ON` and linked against the active Python environment.

Using CMake Presets (Recommended):
```bash
# Configure, build, and test via single workflow preset:
cmake --workflow --preset validation

# Or step-by-step:
cmake --preset validation
cmake --build --preset validation
ctest --preset validation
```

Or using manual flags:
```bash
cmake -B build \
  -DPython3_EXECUTABLE=$(which python3) \
  -DJGAP_BUILD_VALIDATION=ON
cmake --build build -j8
```

---

## 3. Running Validation Tests

### 3.1 Running via CTest
Using the CTest preset:
```bash
ctest --preset validation
```

Or running directly against a custom build directory:
```bash
ctest --test-dir build -L validation --output-on-failure
```

To run a specific validation test:
```bash
ctest --test-dir build -R ValidationElementalProperties --output-on-failure
ctest --test-dir build -R ValidationQrQuick --output-on-failure
ctest --test-dir build -R ValidationQRGapFit --output-on-failure
ctest --test-dir build -R ValidationGapEnergyEval --output-on-failure
ctest --test-dir build -R ValidationTabulation --output-on-failure
ctest --test-dir build -R ValidationTabGapEnergyEval --output-on-failure
```

### 3.2 Running Individual Binaries Directly
Each validation component can also be invoked as an independent command-line binary from the repository root:

```bash
# 1. QR variants benchmark (quick mode)
./build/test/validation/ValidateQrVariants --quick

# 2. Coefficient comparison against reference QUIP XMLs
./build/test/validation/ValidateQRGapFit

# 3. Energy and force predictions vs QUIP predictions
./build/test/validation/ValidateGapEnergyEval

# 4. Tabulation spline tables vs reference tabGAP tables
./build/test/validation/ValidateTabulation

# 5. TabGap potential engine vs reference LAMMPS predictions
./build/test/validation/ValidateTabGapEnergyEval

# 6. Physical material properties (refits via standard_gap_fit + ASE + elastic)
python3 test/validation/test_elemental_properties.py
```

---

## 4. Test Components Detailed Breakdown

### 4.1 `ValidateQrVariants` (QR Decomposition Variants Benchmark)
* **Executable**: `./build/test/validation/ValidateQrVariants`
* **Purpose**: Compares solver variants against the mathematical gold standard (`FullQR`) across increasing 3-body basis sizes ($N_{\text{3b}} = 50, 100, 150, 200$):
  * `FullQR`: In-memory Householder QR of the full matrix $A \in \mathbb{R}^{K \times M}$.
  * `BlockIncrementalQR`: Block-by-block sequential QR factorization.
  * `ElementIncrementalQR`: Incremental species-by-species factorization.
  * Out-of-core streaming QR under constrained RAM budgets (`--approx-ram-limit`).
* **Flags**:
  * `--quick`: Tests only $N_{\text{3b}} = 50$, seed 0 for rapid verification.
  * `--csv <path>`: Exports per-case timing and accuracy metrics.
* **Pass Criteria**:
  * $\text{NRMSE} \le 0.01\%$ vs `FullQR`
  * $\text{Cosine Similarity} \ge 0.999999$

### 4.2 `ValidateQRGapFit` (Fit Coefficients vs Reference QUIP Potential)
* **Executable**: `./build/test/validation/ValidateQRGapFit`
* **Purpose**: Verifies that JGAP's QR linear solver reproduces identical regression coefficients when given the exact reference sparse points from a QUIP-generated `gap.xml`.
* **Workflow**:
  1. Parses QUIP XML using `QuipXmlConverter`.
  2. Extracts reference sparse configurations and regularization sigmas.
  3. Fits the JGAP potential on `feni-train.xyz`.
  4. Compares fitted coefficients $c_j$ with QUIP XML coefficients $c_{\text{ref}}$.
* **Pass Criteria**:
  * Overall $\text{NRMSE} \le 0.01\%$
  * Overall $\text{SigRel} \le 0.05\%$
  * Overall $\text{Cosine Similarity} \ge 0.999999$

### 4.3 `ValidateGapEnergyEval` (Energy & Force Predictions vs Reference QUIP XYZ)
* **Executable**: `./build/test/validation/ValidateGapEnergyEval`
* **Purpose**: Verifies that in-engine Gaussian Process evaluation accurately predicts total energy, Cartesian forces, and virials across test configurations.
* **Workflow**:
  1. Converts reference potential XML via `QuipXmlConverter::transform`.
  2. Evaluates energies and forces over test structures from `feni-test.xyz`.
  3. Compares predictions directly against QUIP extended XYZ outputs (`reference_preds/`).
* **Pass Criteria**:
  * Energy $\text{RMSE} \le 0.01\text{ meV/atom}$
  * Force $\text{RMSE} \le 0.05\text{ meV/\AA}$

### 4.4 `ValidateTabulation` (Spline Tables vs Reference tabGAP HDF5)
* **Executable**: `./build/test/validation/ValidateTabulation`
* **Purpose**: Validates JGAP's tabulation pipeline (`utils::standardTabulation`) against tables generated by the original Python `tabGAP` code.
* **Workflow**:
  1. Converts potential XML to `GapPotential`.
  2. Generates `.tabgap.h5` and `.eam.fs` spline grid files into a temporary directory.
  3. Reads 3-body spline datasets (`E_3B_*`) from both generated and reference HDF5 files.
  4. Reads 2-body pair potential and EAM embedding / electron density tables.
* **Pass Criteria**:
  * 3B spline table $\text{NRMSE} \le 0.01\%$
  * EAM embedding & density spline tables exact match within tolerance.

### 4.5 `ValidateTabGapEnergyEval` (TabGapPotential vs Reference LAMMPS Predictions)
* **Executable**: `./build/test/validation/ValidateTabGapEnergyEval`
* **Purpose**: Verifies JGAP's ultra-fast in-memory spline evaluator (`TabGapPotential`) against LAMMPS (`pair_style hybrid/overlay tabgap eam/fs`).
* **Workflow**:
  1. Loads tabulated `.tabgap.h5` and `.eam.fs` into `TabGapPotential`.
  2. Computes total energies and atomic forces across all 528 test structures of `feni-test.xyz`.
  3. Compares against LAMMPS predictions generated with identical cutoff tables.
* **Pass Criteria**:
  * Energy $\text{RMSE} \le 0.01\text{ meV/atom}$
  * Force $\text{RMSE} \le 0.50\text{ meV/\AA}$ (typical observed: $\sim 0.23\text{ meV/\AA}$)

### 4.6 `test_elemental_properties.py` (DFT Physical Material Properties)
* **Script**: `python3 test/validation/test_elemental_properties.py`
* **Purpose**: End-to-end physical validation: refits potentials directly from raw DFT databases via `jgap.standard_gap_fit` + `jgap.standard_tabulation`, evaluates macroscopic materials properties using ASE and `elastic`, and compares with DFT reference literature values.
* **Evaluated Systems**:
  * **Al** (FCC) from `test/structure-databases/db_Al.xyz`
  * **Cu** (FCC) from `test/structure-databases/db_Cu.xyz`
  * **Ni** (FCC) from `test/structure-databases/db_Ni.xyz`
  * **Fe** (BCC) from `test/structure-databases/db_Fe.xyz`
  * **FeNi** (Alloy) from `test/structure-databases/feni-train.xyz`
* **Computed Properties**:
  * Relaxed lattice constant ($a_0$) via `FrechetCellFilter` and BFGS.
  * Bulk cohesive energy ($E_{\text{coh}} = E_{\text{iso}} - E_{\text{bulk}}/N$).
  * Elastic constants ($C_{11}, C_{12}, C_{44}$) and Bulk Modulus ($B$) via elementary strain deformations.
  * Relaxed vacancy formation energy ($E_{\text{vac}}^f$).
* **Caching**:
  * Fitted potentials are cached in `build/validation_pots/`.
  * Subsequent test runs execute in $\sim 10$ seconds using the cached tables.
  * Pass `--refit` to force refitting all potentials from scratch.
* **Tolerances**:
  * Fails with exit code 1 if any property deviates from DFT beyond allowable thresholds ($|a_0 - a_0^{\text{DFT}}| > 0.05$ Å, $|B - B^{\text{DFT}}| > 40$ GPa, etc.).

### 4.7 `analyze_validation.py` (Reporting & Publication Tables)
* **Script**: `python3 test/validation/analyze_validation.py`
* **Purpose**: Post-processing tool to parse CSV metrics produced by validation runs and format markdown or LaTeX tables.
* **Usage**:
  ```bash
  python3 test/validation/analyze_validation.py --csv <path_to_metrics.csv> [--latex] [--out-dir figures/]
  ```

---

## 5. Statistical Error Metrics Reference

The validation suite uses consistent mathematical formulations defined in [`common/ValidationUtils.hpp`](common/ValidationUtils.hpp):

* **Normalized Root-Mean-Square Error (NRMSE %)**:
  $$\text{NRMSE} = \frac{\sqrt{\frac{1}{N}\sum_{i=1}^N (y_i - \hat{y}_i)^2}}{y_{\max} - y_{\min}} \times 100\%$$
* **Signal-to-Deviation Relative Error (SigRel %)**:
  $$\text{SigRel} = \frac{\sqrt{\frac{1}{N}\sum_{i=1}^N (y_i - \hat{y}_i)^2}}{\text{std}(y_{\text{ref}})} \times 100\%$$
* **Cosine Similarity**:
  $$\text{CosineSim} = \frac{\mathbf{y} \cdot \hat{\mathbf{y}}}{\|\mathbf{y}\|_2 \|\hat{\mathbf{y}}\|_2}$$
* **Root-Mean-Square Error (RMSE)**:
  $$\text{RMSE} = \sqrt{\frac{1}{N}\sum_{i=1}^N (y_i - \hat{y}_i)^2}$$
