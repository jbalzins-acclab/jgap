# JGAP

A high-performance C++23 library, command-line tool, and Python framework for fitting and evaluating **Gaussian Approximation Potentials (GAP)**, **Tabulated GAP (tabGAP)**, and **Embedded Atom Method (EAM)** potentials.

---

## 1. Prerequisites

* **CMake $\ge$ 3.25** and **Ninja** (recommended).
* **C++23 Compliant Compiler**:
  * **GCC $\ge$ 15** or **Clang $\ge$ 19** (recommended for `#embed` support of built-in screening parameters).
  * **AppleClang $\ge$ 16** (Xcode 16+) is fully supported.
  * **GCC 14** / older C++23 compilers work seamlessly by loading runtime screening tables from `resources/`.

---

## 2. Dependencies Overview

`jgap` uses a modern hybrid dependency model:

### Automatic Dependencies (via CMake `FetchContent`)
The following C++ libraries are **automatically downloaded and configured** during CMake build. You do **not** need to install them manually:
* **Eigen3** ($\ge 3.4.0$) — Linear algebra template library.
* **HighFive** ($\ge 3.0.0$) — Header-only modern C++ wrapper for HDF5.
* **pugixml** ($\ge 1.15$) — XML parser for QUIP potential conversion (`jgap_convert`).
* **oneTBB** ($\ge 2022.0.0$) — Auto-fetched for multi-core parallelization if not found on the host system.
* **GoogleTest** ($\ge 1.14$) — Unit testing framework (Debug builds only).

### Host System Dependencies
The following native runtime libraries should be present on your host system:
1. **HDF5** (`libhdf5`) — **Required** for reading/writing `.jgap.h5` and `.tabgap.h5` files.
2. **BLAS / OpenBLAS** — **Strongly Recommended** for accelerated linear algebra (`EIGEN_USE_BLAS`). (On macOS, Apple Accelerate is used automatically if OpenBLAS is not present).
3. **oneTBB** (`libtbb`) — **Optional / Auto-fetched** (`HAS_TBB`). If a host/module TBB is present, it will be used; otherwise CMake automatically downloads and compiles oneTBB.
4. **Python $\ge$ 3.10 + pybind11** — **Optional** for building the `jgap` Python package and ASE calculator.

---

## 3. How to Check Existing Host Dependencies

Before installing new packages, you can verify whether your system or HPC cluster already provides them:

### Using `pkg-config`
```bash
pkg-config --modversion hdf5 openblas tbb
```

### Checking Package Managers
* **macOS (Homebrew)**:
  ```bash
  brew list --formula | grep -E 'hdf5|openblas|tbb'
  ```
* **Debian / Ubuntu**:
  ```bash
  dpkg -l | grep -E 'libhdf5-dev|libopenblas-dev|libtbb-dev'
  ```
* **Conda / Mamba**:
  ```bash
  conda list | grep -E 'hdf5|openblas|tbb'
  ```

### Checking HPC Environment Modules (`Lmod` / `module spider`)
On supercomputing clusters, modules often require loading prerequisite compilers or MPI stacks first. Use `module spider` to inspect available modules and their prerequisite chains:

```bash
# Check compiler, OpenBLAS, and HDF5
module spider gcc
module spider openblas
module spider hdf5

# Check TBB (try common module aliases if 'tbb' is not found)
module spider tbb
module spider onetbb
module spider intel-oneapi-tbb
module spider imkl
```

Example workflow on an Lmod-based HPC cluster:
```bash
# 1. Load compiler and prerequisite stacks (as indicated by module spider)
module load gcc/15.2.0 openmpi/5.0.10

# 2. Load math and I/O libraries
module load openblas/0.3.30 hdf5/1.14.6
```

---

## 4. Installing Host Dependencies

If dependencies are missing on your workstation or cluster, install them using your preferred method:

### macOS (Homebrew)
```bash
brew install cmake ninja hdf5 openblas tbb
```
*(Apple Accelerate is also detected automatically out-of-the-box on macOS).*

### Ubuntu / Debian (`apt`)
```bash
sudo apt update
sudo apt install -y cmake ninja-build build-essential \
                    libhdf5-dev libopenblas-dev libtbb-dev \
                    python3-dev python3-pip
```

### Fedora / RHEL (`dnf`)
```bash
sudo dnf install -y cmake ninja-build gcc-c++ \
                    hdf5-devel openblas-devel tbb-devel \
                    python3-devel
```

### Arch Linux (`pacman`)
```bash
sudo pacman -S cmake ninja hdf5 openblas intel-oneapi-tbb python
```

### Conda / Mamba (Recommended for User-Space HPC Environments)
```bash
conda install -c conda-forge cmake ninja compilers \
                            hdf5 openblas tbb-devel pybind11
```
*(When using Conda, CMake will automatically locate dependencies inside `$CONDA_PREFIX`).*

> **Note on oneTBB**: If your environment does not have TBB installed, CMake will automatically fetch and compile `oneTBB` via `FetchContent` during configuration. No manual installation is needed!

---

## 5. Building and Testing with CMake Presets

`jgap` provides built-in `CMakePresets.json` profiles for all standard workflows:

### Fast Developer Workflow (Debug + Build + Run All Tests)
```bash
cmake --workflow --preset dev
```

### Release Build (Optimized with `-O3 -ffast-math -march=native`)
```bash
cmake --preset release
cmake --build --preset release
```

### Install Library, CLI & Python Bindings (to `$HOME/.local`)
```bash
cmake --workflow --preset install
```
*(To install to a custom prefix or virtual environment, run `cmake --preset release -DCMAKE_INSTALL_PREFIX=/path/to/prefix && cmake --build --preset install`).*

### Run Unit Tests
```bash
ctest --preset debug
```

### AddressSanitizer (Memory Diagnostics)
```bash
cmake --preset asan
cmake --build --preset asan
ctest --preset asan
```

---

## 6. Python Package & ASE Calculator

After building or installing `jgap`, the Python bindings are available in `python/jgap`:

```python
import jgap
from ase.io import read

# Standard GAP potential fit
fitter = jgap.StandardGapFit(
    cutoff_2b=5.0,
    cutoff_3b=4.0,
    delta_2b=0.01,
    delta_3b=0.05,
    approx_ram_limit_gb=16.0
)
fitter.fit("train.xyz")
fitter.save("potential.jgap.h5")

# Tabulate into fast EAM / tabGAP
jgap.StandardTabulation.tabulate("potential.jgap.h5", "potential")
```

Using as an ASE Calculator:
```python
from jgap.ase import JGAPCalculator
from ase.io import read

atoms = read("structure.xyz")
atoms.calc = JGAPCalculator("potential.jgap.h5")

energy = atoms.get_potential_energy()
forces = atoms.get_forces()
```

---

## 7. Command-Line Tools

* **Predict Energy & Forces**:
  ```bash
  jgap --predict potential.jgap.h5 input.xyz output.xyz
  ```
* **Tabulate Potential (to EAM `.eam.fs` and tabGAP `.tabgap.h5`)**:
  ```bash
  jgap --tabulate potential.jgap.h5
  ```
* **Convert QUIP XML to HDF5**:
  ```bash
  jgap_convert potential.xml potential.jgap.h5
  ```