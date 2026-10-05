# Testing in JGAP

`jgap` utilizes a two-tier testing architecture designed to balance ultra-fast local feedback during routine development with rigorous, full-scale scientific validation against external simulation engines and experimental/DFT reference baselines.

---

## 1. Testing Architecture Overview

```
test/
├── unit/                 # Tier 1: Fast GoogleTest unit test suites (jgap_tests)
│   ├── core/             # Atoms, Matrix, cutoffs, transformations, sparsification, splines
│   ├── impl/fit/         # QRGapFit, BlockIncrementalQR, ElementIncrementalQR
│   └── serialization/    # HDF5 serialization round-trip verification
├── validation/           # Tier 2: Automated end-to-end scientific validation suite
│   ├── ValidateQrVariants.cpp       # QR solver variants benchmark
│   ├── ValidateQRGapFit.cpp         # QUIP XML coefficient reproduction
│   ├── ValidateGapEnergyEval.cpp    # In-engine energy & force prediction verification
│   ├── ValidateTabulation.cpp       # Spline table generation vs tabGAP HDF5
│   ├── ValidateTabGapEnergyEval.cpp # TabGapPotential engine vs LAMMPS
│   ├── test_elemental_properties.py # End-to-end DFT physical property validation
│   └── ReadMe.md                    # Detailed documentation for validation suite
└── resources/            # Test datasets, sample structures, and reference potentials
    ├── structure-databases/ # Training and test databases (feni-train.xyz, feni-test.xyz, db_*.xyz)
    ├── reference/        # Non-zip reference potentials and benchmark predictions for unit tests
    ├── quip_potentials/  # Baseline XML potentials
    ├── validation.zip    # Compressed reference validation archives (52 MB)
    └── validation/       # Extracted reference baselines (QUIP, tabGAP, LAMMPS, DFT)
```

---

## 2. Tier 1: Fast Unit Tests (`jgap_tests`)

### Purpose
Unit tests verify the correctness, numerical stability, and edge cases of individual components, mathematical primitives, and data structures:
* **Atomic Data Structures**: `Atoms`, `Matrix`, neighbor list construction, periodic boundary wrapping.
* **Transformations & Kernels**: 2-body, 3-body, and EAM transformations, squared exponential kernels, permutation invariance.
* **Sparsification**: `HistogramUniformSparsifier` distribution and binning.
* **Linear Solvers & Fitters**: `QRGapFit`, `BlockIncrementalQRGapFit`, and `ElementIncrementalQRGapFit` mathematical equivalence and low-memory streaming.
* **Splines & Tabulation**: 1D Natural, Hermite, and B-Spline interpolation, 3D cubic B-splines, and `TabGapPotential` evaluation.
* **Serialization**: Complete HDF5 serialization and deserialization round-trips for all potentials, components, and transformations.

### Building & Running Unit Tests
Unit tests are automatically discovered from all `test/unit/**/*.cpp` files and compiled into the `jgap_tests` executable using GoogleTest:

```bash
# Build the unit tests target
cmake --build build --target jgap_tests -j8

# Run all unit tests directly
./build/jgap_tests

# Or run via CTest filtering
ctest --test-dir build -R jgap_tests --output-on-failure

# Run a specific test suite or test filter
./build/jgap_tests --gtest_filter="QrGapFits.*"
./build/jgap_tests --gtest_filter="TestTabGapPotential.*"
```

---

## 3. Tier 2: Automated Validation Suite (`ctest -L validation`)

### Purpose
The validation suite conducts end-to-end verification against external reference software without requiring those external binaries to be installed at test time:
1. **`ValidateQrVariants`**: Benchmarks full QR, block incremental QR, and element incremental QR across scaling basis dimensions ($M$).
2. **`ValidateQRGapFit`**: Verifies exact coefficient reproduction ($c_j$) against reference QUIP XML potentials.
3. **`ValidateGapEnergyEval`**: Validates total energy, Cartesian forces, and virials against reference QUIP predictions.
4. **`ValidateTabulation`**: Compares JGAP spline tabulation grid values against reference Python `tabGAP` HDF5 tables.
5. **`ValidateTabGapEnergyEval`**: Validates `TabGapPotential` energy and force evaluation against LAMMPS predictions.
6. **`ValidationElementalProperties`**: Refits potentials directly from raw DFT databases via `jgap.standard_gap_fit`, evaluates material properties ($a_0, B, C_{ij}, E_{\text{vac}}^f, E_{\text{coh}}$) via ASE + `elastic`, and compares with DFT literature references.

> [!TIP]
> For in-depth workflows, input files, pass/fail thresholds, and analysis tools for Tier 2, refer to the dedicated **[Validation Suite ReadMe](validation/ReadMe.md)**.

### Running Validation Benchmarks
First, extract the reference validation data (archived in `test/resources/validation.zip`):
```bash
unzip -q test/resources/validation.zip -d test/resources/
```

Run validation tests via CMake Presets:
```bash
# Full workflow (configure, build, test):
cmake --workflow --preset validation

# Or run tests using CTest label:
ctest --test-dir build -L validation --output-on-failure
```

---

## 4. Test Resources & Reference Data

* **`test/structure-databases/`**: Contains raw training databases (`feni-train.xyz`, `feni-test.xyz`, `db_Al.xyz`, `db_Cu.xyz`, `db_Ni.xyz`, `db_Fe.xyz`).
* **`test/resources/reference/`**: Non-zip reference baselines tracked in git for fast unit test execution (`reference_pots/`, `reference_preds/`, `reference_tables/`, `reference_lammps/`).
* **`test/resources/validation.zip`**: Contains precomputed reference baselines:
  * `reference_pots/`: Reference QUIP XML potentials (`gap.xml`).
  * `reference_preds/`: QUIP extended XYZ predictions with reference energies and forces.
  * `reference_tables/`: Python `tabGAP` HDF5 spline tables (`.tabgap.h5`, `.eam.fs`).
  * `reference_lammps/`: Reference LAMMPS predictions.
  * `elemental_property_comparison.csv`: Reference DFT lattice and elastic constants.
