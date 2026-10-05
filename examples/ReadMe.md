# JGAP Examples & Standalone Usage

This directory contains standalone examples demonstrating how to use `jgap` from both C++ and Python for potential fitting, tabulation, and simulation with ASE. All fitting examples use `ElementIncrementalQRGapFit` with a default 2 GB RAM limit (adjustable via `--ram-limit`).

---

## 1. Directory Overview

| Directory | Language | Description |
| :--- | :--- | :--- |
| **[`basic_fit/`](basic_fit/BasicFit.cpp)** | C++ | Manual component assembly without helper utilities: 2b + 3b + EAM GAP with `SquaredExpKernel`, fast $12 \times 12 \times 12$ 3B tabulation grid, and `ElementIncrementalQRGapFit`. |
| **[`custom_fit/`](custom_fit/CustomFit.cpp)** | C++ | Advanced multi-component Fe-Ni potential showcasing multiple kernel types (`WendlandKernel`, `CauchyKernel`, `SquaredExpKernel`), MEAM descriptor (`ThreeBodySum` + `MeamTransformation`), `PerriotPolynomialCutoff`, and `ScaledRegularizationRules`. |
| **[`custom_experiment_fit/`](custom_experiment_fit/CustomExperimentFit.cpp)** | C++ | Demonstrates user-defined custom kernels: implements `FractionalExpKernel` ($-\|r/\ell\|^{1.5}$, strictly positive definite by Schoenberg's theorem), 2b-only GAP, directly tabulates to `.tabgap.h5` without saving GAP potential. |
| **[`standard_fit/`](standard_fit/)** | C++ & Python | Standard fitting and tabulation pipeline with `ElementIncrementalQRGapFit`: [C++ version (`StandardFit.cpp`)](standard_fit/StandardFit.cpp) and [Python version (`standard_fit.py`)](standard_fit/standard_fit.py). |
| **[`ase_integration/`](ase_integration/test_potetnial_with_ase.py)** | Python | Evaluates fitted potentials with ASE: per-config-type validation errors (Energy RMSE/MAE in meV/atom, Force RMSE in meV/Å, Virial RMSE in meV/atom) and bulk Fe properties with capped relaxation steps (`--max-steps 100`). |

---

## 2. Compiling Standalone C++ Examples (Without CMake)

Assuming `jgap` has been installed (e.g. via `cmake --workflow --preset install`), you can compile standalone C++ programs directly with `c++`:

### Direct Compiler Command Line

```bash
# Basic fit (manual component assembly)
c++ -std=c++23 -O3 -march=native examples/basic_fit/BasicFit.cpp -ljgap -o examples/basic_fit/basic_fit

# Standard fit
c++ -std=c++23 -O3 -march=native examples/standard_fit/StandardFit.cpp -ljgap -o examples/standard_fit/standard_fit

# Custom fit (multi-kernel, MEAM descriptor, Perriot cutoff)
c++ -std=c++23 -O3 -march=native examples/custom_fit/CustomFit.cpp -ljgap -o examples/custom_fit/custom_fit

# Custom experiment fit (user-defined FractionalExpKernel, direct tabulation)
c++ -std=c++23 -O3 -march=native examples/custom_experiment_fit/CustomExperimentFit.cpp -ljgap -o examples/custom_experiment_fit/custom_experiment_fit
```

> [!NOTE]
> If installed to a non-standard prefix (such as `$HOME/.local`), ensure your compiler search paths include the include and library directories or pass `-I$HOME/.local/include -L$HOME/.local/lib -Wl,-rpath,$HOME/.local/lib`.

---

## 3. Running the Examples

### Running C++ Fit Executables

All fitting examples default to a 2 GB RAM limit, which can be adjusted with `--ram-limit <gb>`:

```bash
# Basic fit with 2 GB RAM limit (default)
./examples/basic_fit/basic_fit test/resources/structure-databases/db_Fe.xyz fe_basic

# Custom RAM limit (e.g. 4.0 GB)
./examples/basic_fit/basic_fit test/resources/structure-databases/db_Fe.xyz fe_basic --ram-limit 4.0

# Standard fit
./examples/standard_fit/standard_fit test/resources/structure-databases/db_Fe.xyz fe_standard --ram-limit 2.0

# Custom multi-kernel Fe-Ni fit with MEAM descriptor
./examples/custom_fit/custom_fit test/resources/structure-databases/feni-train.xyz feni-custom --ram-limit 2.0

# Custom experiment fit (user-defined kernel, direct tabulation to .tabgap.h5)
./examples/custom_experiment_fit/custom_experiment_fit test/resources/structure-databases/db_Fe.xyz fe_exp --ram-limit 2.0
```

### Running Python Examples

Ensure your environment with `jgap` is active:

```bash
# Run standard fit on sample dataset (defaults to 2 GB RAM limit)
python examples/standard_fit/standard_fit.py \
    test/resources/structure-databases/db_Fe.xyz \
    fe_potential \
    --ram-limit 2.0

# Evaluate potential with ASE (per-config-type E/F/Virial RMSE and bulk BCC Fe properties)
python examples/ase_integration/test_potetnial_with_ase.py \
    fe_exp.tabgap.h5 \
    test/resources/structure-databases/db_Fe.xyz

# Optionally disable per-config-type error table:
python examples/ase_integration/test_potetnial_with_ase.py \
    fe_exp.tabgap.h5 \
    test/resources/structure-databases/db_Fe.xyz \
    --no-config-errors
```
