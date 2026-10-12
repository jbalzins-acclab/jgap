# JGAP Examples & Standalone Usage

This directory contains standalone examples demonstrating how to use `jgap` from both C++ and Python for potential fitting, tabulation, and simulation with ASE.

---

## 1. Directory Overview

| Directory | Language | Description |
| :--- | :--- | :--- |
| **[`basic_fit/`](basic_fit/BasicFit.cpp)** | C++ | Manual 2b + 3b + EAM GAP component assembly without helper utilities. |
| **[`custom_fit/`](custom_fit/CustomFit.cpp)** | C++ | Advanced multi-component potential fitting, showcasing multiple kernel types (`WendlandKernel`, `CauchyKernel`, `SquaredExpKernel`), MEAM descriptor (`ThreeBodySum` + `MeamTransformation`), `PerriotPolynomialCutoff`, and `ScaledRegularizationRules`. |
| **[`custom_experiment_fit/`](custom_experiment_fit/CustomExperimentFit.cpp)** | C++ | Demonstrates how to use user-define a custom kernel for fast testing. |
| **[`standard_fit/`](standard_fit/)** | C++ & Python | Standard fitting and tabulation pipeline: [C++ version (`StandardFit.cpp`)](standard_fit/StandardFit.cpp) and [Python version (`standard_fit.py`)](standard_fit/standard_fit.py). |
| **[`ase_integration/`](ase_integration/test_potential_with_ase.py)** | Python | Evaluates fitted potentials with ASE: bulk properties and test-set errors. |

---

## 2. Compiling Standalone C++ Examples (Without CMake)

Assuming `jgap` has been installed (e.g. via `cmake --workflow --preset install`), you can compile standalone C++ programs directly with `c++`(=`g++` or `clang++`):

### Direct Compiler Command Line

```bash
c++ -std=c++23 -O3 -march=native examples/basic_fit/BasicFit.cpp -ljgap -o examples/basic_fit/basic_fit
c++ -std=c++23 -O3 -march=native examples/standard_fit/StandardFit.cpp -ljgap -o examples/standard_fit/standard_fit
c++ -std=c++23 -O3 -march=native examples/custom_fit/CustomFit.cpp -ljgap -o examples/custom_fit/custom_fit
c++ -std=c++23 -O3 -march=native examples/custom_experiment_fit/CustomExperimentFit.cpp -ljgap -o examples/custom_experiment_fit/custom_experiment_fit
```

> [!NOTE]
> If installed to a non-standard prefix (such as `$HOME/.local` or virtaul environment's prefix), ensure your compiler search paths include the include and library directories or pass `-I$HOME/.local/include -L$HOME/.local/lib -Wl,-rpath,$HOME/.local/lib`.

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

# Evaluate bulk properties only (testset.xyz is optional):
python examples/ase_integration/test_potential_with_ase.py fe_exp.tabgap.h5

# Evaluate bulk properties and per-config-type E/F/Virial errors:
python examples/ase_integration/test_potential_with_ase.py \
    fe_exp.tabgap.h5 \
    test/resources/structure-databases/db_Fe.xyz
```
