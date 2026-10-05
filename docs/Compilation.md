# Compilation and Installation Guide

This guide details how to build, install, and configure `jgap` as both a high-performance C++23 shared library and a Python package integrated into your virtual environment of choice.

---

## 1. Prerequisites

### C++23 Compiler
`jgap` utilizes modern C++23 language features (including concepts, standard library extensions, and `#embed` where supported).
* **GCC $\ge$ 15** or **Clang $\ge$ 19** (recommended; provides `#embed` support for built-in Coulomb screening datasets).
* **AppleClang $\ge$ 16** (Xcode 16+) is fully supported on macOS.
* **GCC 14** works seamlessly by falling back to loading screening tables from `resources/` at runtime.

### Build Tools
* **CMake $\ge$ 3.25**
* **Ninja** (recommended generator for fast parallel builds)

---

## 2. Dependencies Overview

### Automatically Downloaded Dependencies (CMake `FetchContent`)
The following lightweight dependencies are downloaded and configured automatically by CMake during the build process:
* **Eigen3** ($\ge 3.4.0$) — High-performance template library for linear algebra.
* **HighFive** ($\ge 3.0.0$) — Modern C++ header-only interface to HDF5.
* **pugixml** ($\ge 1.15$) — Fast XML parser for QUIP potential conversion.
* **GoogleTest** ($\ge 1.14$) — Unit testing framework (built only when tests are enabled).

### Host System Dependencies
Ensure the following native libraries are installed on your host system:
1. **HDF5** (`libhdf5`) — **Required** for reading/writing `.jgap.h5` and `.tabgap.h5` potential files.
2. **BLAS / OpenBLAS** — **Strongly Recommended** for accelerating linear algebra operations (`EIGEN_USE_BLAS`). On macOS, Apple Accelerate is automatically detected if OpenBLAS is not present.
3. **OpenMP** — **Recommended** for multi-threaded energy evaluation, force evaluation, and neighbor list construction (`HAS_OPENMP`). Built into GCC/Clang/Intel; on macOS available via Homebrew (`libomp`).
4. **Python $\ge$ 3.10 & pybind11** — **Required for Python bindings** (`jgap` Python package and ASE calculator).

### Installing Host Dependencies by Operating System

#### macOS (Homebrew)
```bash
brew install cmake ninja hdf5 openblas libomp
```
*(Apple Accelerate is also detected automatically out-of-the-box on macOS).*

#### Ubuntu / Debian (`apt`)
```bash
sudo apt update
sudo apt install -y cmake ninja-build build-essential \
                    libhdf5-dev libopenblas-dev libomp-dev \
                    python3-dev python3-pip
```

#### Fedora / RHEL (`dnf`)
```bash
sudo dnf install -y cmake ninja-build gcc-c++ \
                    hdf5-devel openblas-devel libgomp \
                    python3-devel
```

#### Arch Linux (`pacman`)
```bash
sudo pacman -S cmake ninja hdf5 openblas openmp python
```

#### Conda / Mamba (Recommended for User-Space HPC Environments)
```bash
conda install -c conda-forge cmake ninja compilers \
                            hdf5 openblas pybind11
```

---

## 3. Targeting Your Python Virtual Environment of Choice

To install `jgap` directly into your active Python environment (whether created via `venv`, `conda`, `poetry`, or `uv`), explicitly supply the Python executable and pybind11 directory to CMake:

### Step 1: Activate Your Environment
```bash
# Example with venv
source ~/my_env/bin/activate

# Or example with Conda
conda activate my_env
```

### Step 2: Ensure `pybind11` is Installed in the Environment
```bash
pip install pybind11 numpy ase
```

### Step 3: Configure and Install
Point CMake to your environment's Python interpreter:
```bash
# Configure Release build targeting the active virtual environment
cmake --preset release \
  -DPython3_EXECUTABLE=$(which python) \
  -DCMAKE_INSTALL_PREFIX=$VIRTUAL_ENV

# Build and install C++ library, headers, CLI, and Python package
cmake --build --preset install
```

> [!NOTE]
> If you are using Conda, replace `-DCMAKE_INSTALL_PREFIX=$VIRTUAL_ENV` with `-DCMAKE_INSTALL_PREFIX=$CONDA_PREFIX`. If installing to your user space outside a venv, default to `-DCMAKE_INSTALL_PREFIX=$HOME/.local`.

### Verifying the Python Installation
```bash
python -c "import jgap; print(f'jgap loaded from {jgap.__file__}')"
```

### In-Tree Development (Without Installing)
When CMake builds the Python extension, the compiled `_jgap` shared module is automatically copied directly into `python/jgap/`. You can use it immediately without installing by setting `PYTHONPATH`:
```bash
export PYTHONPATH="$(pwd)/python:$PYTHONPATH"
python -c "import jgap; print(jgap)"
```

---

## 4. CMake Presets Reference

`jgap` provides declarative presets in `CMakePresets.json`:

| Preset Name | Type | Description |
| :--- | :--- | :--- |
| `debug` | Configure | Debug build (`-g`), unit tests enabled, compile commands exported. |
| `release` | Configure | Optimized release build (`-O3 -ffast-math -march=native`). |
| `release-tests` | Configure | Optimized release build with unit tests enabled (`JGAP_BUILD_TESTS=ON`). |
| `relwithdebinfo` | Configure | Optimized release with debug symbols (`-O3 -g -march=native`). |
| `asan` | Configure | Debug build with AddressSanitizer & UndefinedBehaviorSanitizer. |
| `validation` | Configure | Optimized release build with Tier 2 validation suite enabled. |
| `dev` | Workflow | Full developer loop: configure `debug`, build `debug`, run unit tests. |
| `ci` | Workflow | CI pipeline: configure `release-tests`, build `release-tests`, run tests. |
| `install` | Workflow | Release workflow: configure `release` and install target. |
| `validation` | Workflow | Configure, build, and run the automated validation benchmarks. |

### Common Preset Commands
```bash
# Fast developer build & unit test run:
cmake --workflow --preset dev

# Release build:
cmake --preset release
cmake --build --preset release

# Install to default prefix ($HOME/.local):
cmake --workflow --preset install

# Run AddressSanitizer checks:
cmake --preset asan
cmake --build --preset asan
ctest --preset asan
```

---

## 5. HPC Cluster Deployment (`Lmod` / Modules)

On supercomputers and clusters, native compilers and optimized math libraries are typically provided through modules.

### Example Workflow
```bash
# 1. Load GCC 15+ or Clang 19+ and OpenBLAS / HDF5 modules
module load gcc/15.2.0
module load openblas/0.3.30
module load hdf5/1.14.6

# 2. Activate your target Python environment
module load python/3.13
source ~/venvs/jgap-env/bin/activate

# 3. Configure and compile
cmake -B build -G Ninja \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_CXX_FLAGS="-O3 -ffast-math -march=native" \
  -DPython3_EXECUTABLE=$(which python) \
  -DCMAKE_INSTALL_PREFIX=$VIRTUAL_ENV

cmake --build build -j16
cmake --install build
```

---

## 6. Troubleshooting & FAQs

### Python Extension Links to Wrong Python Version
If CMake selects a system Python (e.g., Python 3.14) rather than your virtual environment's Python (e.g., Python 3.13):
1. Delete the `build/` directory: `rm -rf build`.
2. Explicitly specify the interpreter path:
   ```bash
   cmake --preset release -DPython3_EXECUTABLE=$VIRTUAL_ENV/bin/python
   ```

### Shared Library Search Paths (`RPATH`)
If running a C++ executable returns `error while loading shared libraries: libjgap.so`:
* On Linux, ensure your install directory is in `LD_LIBRARY_PATH`:
  ```bash
  export LD_LIBRARY_PATH=$CMAKE_INSTALL_PREFIX/lib:$LD_LIBRARY_PATH
  ```
* On macOS, add to `DYLD_LIBRARY_PATH`:
  ```bash
  export DYLD_LIBRARY_PATH=$CMAKE_INSTALL_PREFIX/lib:$DYLD_LIBRARY_PATH
  ```
