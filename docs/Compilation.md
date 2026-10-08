# Compilation & Installation Guide

`jgap` is a high-performance C++23 library, CLI toolkit, and Python framework. Because different environments have different privileges, package managers, and hardware capabilities, dedicated installation guides and workflow presets are detailed below.

---

## 1. Quickstart: One-Command Build & Install

If your system already has the prerequisite C++23 compiler and libraries installed, you can configure, build, and install the library, CLI tools, and Python interface in a single command using CMake's unified workflow preset:

```bash
# Configure, build optimized Release targets, and install in a single command
cmake --workflow --preset install
```

> [!TIP]
> **Active Python Virtual Environment**: To install into your virtualenv or Conda environment, pass `-DCMAKE_INSTALL_PREFIX=$VIRTUAL_ENV` (or `$CONDA_PREFIX`):
> ```bash
> cmake --preset release -DCMAKE_INSTALL_PREFIX=$VIRTUAL_ENV -DPython3_EXECUTABLE=$(which python)
> cmake --build --preset install
> ```

---

## 2. Quick Selector: Choose Your Environment

Select the detailed guide matching your platform and access level:

| Platform | Access Level / Tool | Guide |
| :--- | :--- | :--- |
| **macOS** (Apple Silicon) | Homebrew available | [**macOS Installation Guide**](installation/macos.md) |
| **Linux** (Ubuntu, Debian, Fedora, Arch) | `sudo` root rights available | [**Linux (`sudo`) Guide**](installation/linux_sudo.md) |
| **Linux** (Any distribution / Server) | No `sudo` — Conda / Mamba / Pixi | [**Linux (Conda / User-Space) Guide**](installation/linux_conda.md) |
| **Linux** (Any distribution / HPC) | No `sudo` — Apptainer / Singularity `.sif` | [**Linux (Apptainer `.sif`) Guide**](installation/linux_apptainer.md) |

---

## 3. Core Technical Requirements Summary

* **Compilers**:
  * **Linux**: **GCC $\ge$ 15** or **Clang $\ge$ 19** (GCC 14 does not support `#embed`).
  * **macOS**: AppleClang $\ge$ 16 (Xcode 16+ on Apple Silicon) or Homebrew LLVM Clang.
* **Build System**: CMake $\ge$ 3.25 and Ninja.
* **C++ Libraries**: HDF5 (`libhdf5`), OpenBLAS (or Apple Accelerate), OpenMP.
* **Python**: Python $\ge$ 3.9 with `pybind11`, `numpy`, and `ase`.

---

## 4. HPC Clusters & `-march=native` Compilation

Pre-compiled binary wheels distributed on PyPI target generic CPU architectures for universal portability. However, in High-Performance Computing (HPC) workflows—where descriptor evaluation, neighbor listing, and 3D B-spline tabulation run millions of times—hardware-level CPU vectorization is critical.

Compiling `jgap` directly on your HPC cluster compute nodes enables **`-march=native`** vector extensions (e.g. **AVX-512**, **AVX2**, and **FMA**), yielding substantial performance gains.

### HPC Step-by-Step Instructions

1. **Load Environment Modules** (example for SLURM / Lmod clusters with GCC 15):
   ```bash
   module purge
   module load gcc/15     # Or llvm/clang-19
   module load openblas
   module load hdf5
   module load cmake
   module load ninja
   module load python/3.12
   ```

2. **Activate Your Python Environment**:
   ```bash
   python3 -m venv ~/jgap-env
   source ~/jgap-env/bin/activate
   pip install --upgrade pip
   pip install pybind11 numpy scipy ase
   ```

3. **Compile and Install with `-march=native`**:
   The `release` preset applies `-O3 -ffast-math -march=native` by default:
   ```bash
   # Configure directly on a compute node (or interactive job) to target its exact microarchitecture:
   cmake --preset release \
     -DCMAKE_INSTALL_PREFIX=$VIRTUAL_ENV \
     -DPython3_EXECUTABLE=$(which python)

   # Build in parallel and install
   cmake --build --preset install -j$(nproc)
   ```

4. **Verify SIMD Vectorization in SLURM Job Scripts**:
   ```bash
   #!/bin/bash
   #SBATCH --job-name=jgap_fit
   #SBATCH --nodes=1
   #SBATCH --cpus-per-task=16
   #SBATCH --time=04:00:00

   source ~/jgap-env/bin/activate
   export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

   jgap --predict potential.jgap.h5 structures.xyz output.xyz
   ```

---

## 5. CMake Workflow Presets Reference

`jgap` provides unified CMake presets (`CMakePresets.json`) across all platforms:

| Preset Name | Type | Description |
| :--- | :--- | :--- |
| `release` | Configure | Optimized release configuration (`-O3 -ffast-math -march=native`). |
| `release-tests` | Configure | Release build with unit tests enabled (`JGAP_BUILD_TESTS=ON`). |
| `debug` | Configure | Debug symbols, compile commands export (`compile_commands.json`). |
| `validation` | Configure | Release configuration with automated validation test suite. |
| `asan` | Configure | AddressSanitizer and UndefinedBehaviorSanitizer instrumentation. |
| `install` | Build | Builds Release targets and installs them to `CMAKE_INSTALL_PREFIX`. |

### Execution Commands
```bash
# Build and install Release binaries in a single step (Recommended)
cmake --workflow --preset install

# Run unit tests in Release mode
cmake --workflow --preset ci

# Build and execute validation test suite
cmake --workflow --preset validation

# Debug build and test execution
cmake --workflow --preset dev

# AddressSanitizer run
cmake --preset asan && cmake --build --preset asan && ctest --preset asan
```

---

## 6. Troubleshooting & FAQs

### Python Extension Links to Wrong Python Version
If CMake selects a system Python rather than your active virtual environment's Python:
1. Delete the `build/` directory: `rm -rf build`.
2. Explicitly specify the interpreter path during configuration:
   ```bash
   cmake --preset release -DPython3_EXECUTABLE=$(which python)
   ```

### Pybind11 Not Found During CMake Configuration
If you see `-- pybind11 not found: skipping _jgap Python module build`:
Ensure `pybind11` is installed inside your active environment:
```bash
pip install pybind11
```
and re-run CMake configuration.
